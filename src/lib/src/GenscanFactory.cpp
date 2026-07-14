//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "GenscanFactory.h"
#include "GenePredDataSet.h"
#include "parse_help.h"
#include "SpliceSiteDataSet.h"
#include "genemodels.h"
#include "DnaStr.h"
#include <fstream>
#include <stdexcept>

using std::endl;
using std::cerr;
using std::invalid_argument;
using std::exception;
using std::cout;

static string conv2exon(const string& str) {
  if(strstr(str.c_str(),"Init")) return "initial"; 
  else if(strstr(str.c_str(),"Intr")) return "internal"; 
  else if(strstr(str.c_str(),"Term")) return "terminal"; 
  else if(strstr(str.c_str(),"Sngl")) return "single"; 
  else return "";
}

static void stripExonId(string& str) {
  for(unsigned i = 0; i < str.length(); ++i) {
    if(str[i] == '.') {
      str.erase(i,str.length());
    }
  }
}

DataSetAbstractProduct* GenscanFactory::createDataSetAbstractProduct(const string& filename, const string& type, const string& evidence_type, const string& seq) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw invalid_argument(filename + " is not a valid input file ");
  const int buff_size = 4096;
  char buff[buff_size];
  SeqLoc sl(dsu::eEither,-1,-1);
  // create a new list of genes
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,type);
  int debug_cnt = 0, geneId = 0;
  string geneIdStr,lastGeneIdStr;
  while( ifs.getline(buff,buff_size) ) {
    if( strstr(buff,"Predicted genes/exons:") ) {
      break;
    }
  }
  ifs.getline(buff,buff_size);
  ifs.getline(buff,buff_size);
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '\n' ) {
      continue;
    }
    if( strstr(buff,"Predicted protein") ) {
      break;
    }
    ++debug_cnt;
    string gffid = type;
    string exon_type = get_elem(buff,2);
    exon_type = conv2exon(exon_type);
    if( exon_type == "") {
      continue;
    }
    int end5 = atoi(get_elem(buff,4).c_str());
    int end3 = atoi(get_elem(buff,5).c_str());
    bool validEnd5 = (end5 >= 1 && end5 < (signed)seq.size());
    bool validEnd3 = (end3 >= 1 && end3 < (signed)seq.size());
    if( !validEnd5 || !validEnd3 ) {
      cout<<"Warning: invalid sequence coordinates (ignore line), problem parsing line ["<<buff<<"]";
      if( !validEnd5 ) {
	cout<<" 5' coordinate = "<<end5;
      }
      if( !validEnd3 ) {
	cout<<" 3' coordinate = "<<end3;
      }
      cout<<endl;
      continue;
    }
    double score = atof(get_elem(buff,12).c_str());
    string geneIdStr = get_elem(buff,1);
    stripExonId(geneIdStr);
    if( geneIdStr != lastGeneIdStr ) {
      ++geneId;
    }
    string strnd = get_elem(buff,3);
    dsu::Strand_t parseStrnd = dsu::eEither; 
    if( strnd == "-" ) {
      parseStrnd = dsu::eNeg;
    } else {
      parseStrnd = dsu::ePos;
    }
    try {
      assert(end5 >= 0 && end3 >= 0);
      SeqLoc theSl(parseStrnd,end5,end3);
      newGenePredDataSet->insert(gffid,geneId,theSl,exon_type,score);
    } catch(exception& e) {
      cerr<<e.what()<<endl;
      cerr<<"attempting to continue: ["<<filename<<"]"<<endl;
    }
    lastGeneIdStr = geneIdStr;
  }
  if(!debug_cnt) {
    throw invalid_argument(filename + " : unable to read file");
  }
  return newGenePredDataSet;
}
