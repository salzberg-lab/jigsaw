//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "SnapFactory.h"
#include "GenePredDataSet.h"
#include "SpliceSiteDataSet.h"
#include "genemodels.h"
#include "DnaStr.h"
#include <fstream>
#include <stdexcept>
#include <assert.h>

using std::cout;
using std::endl;

static string conv2exon(const string& str) {
  if(strstr(str.c_str(),"final")) return "terminal"; 
  else if(strstr(str.c_str(),"Eterm")) return "terminal"; 
  else if(strstr(str.c_str(),"term")) return "terminal"; 
  else if(strstr(str.c_str(),"Einit")) return "initial"; 
  else if(strstr(str.c_str(),"initial")) return "initial"; 
  else if(strstr(str.c_str(),"Exon")) return "internal"; 
  else if(strstr(str.c_str(),"internal")) return "internal"; 
  else if(strstr(str.c_str(),"Esngl")) return "single"; 
  else if(strstr(str.c_str(),"single")) return "single"; 
  else if(strstr(str.c_str(),"first")) return "initial"; 
  else if(strstr(str.c_str(),"last")) return "initial"; 
  else if(strstr(str.c_str(),"start")) return "initial"; 
  else if(strstr(str.c_str(),"stop")) return "terminal"; 
  else return "internal";
}

#if 0
static bool useCurationFormat(char* buff) {
  int cnt = 0;
  for(char* ptr = buff; *ptr && *ptr != '\n'; ++ptr) {
    if(*ptr == ' ') ++cnt;
  }
  if(cnt==2) return true; 
  return false;
}
#endif

static string get_tab_elem(const string& str, unsigned ith_tab) {
  string::size_type pos=0,stop=string::npos; 
  assert(pos < str.size());
  for(unsigned i = 0; i <= ith_tab; ++i) {
    stop = str.find("\t",pos);
    if(stop==string::npos)
      stop = str.size();
    if( i < ith_tab ) 
      pos = stop+1;
  }
  string sstr = str.substr(pos,stop-pos);
  return sstr;
}

DataSetAbstractProduct* 
SnapFactory::createDataSetAbstractProduct(const string& filename, const string& type, const string& evidence_type, const string& seq) {
  std::ifstream ifs(filename.c_str());
  if(!ifs) {
    throw std::invalid_argument(filename + " is not a valid input file ");
  }
  const int buff_size = 4096;
  char buff[buff_size];
  SeqLoc sl(dsu::eEither,-1,-1);
  // create a new list of genes
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,type);
  int debug_cnt = 0, geneId = 0;
  string geneIdStr,lastGeneIdStr;
  string gffid = type;
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '#' || *buff == '\n' || *buff == '>' ) 
      continue;
    ++debug_cnt;
    string exon_type = get_tab_elem(buff,0);
    exon_type = conv2exon(exon_type);
    int end5 = atoi(get_tab_elem(buff,1).c_str());
    int end3 = atoi(get_tab_elem(buff,2).c_str());
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
    double score = atof(get_tab_elem(buff,4).c_str());
    string geneIdStr = get_tab_elem(buff,8);
    if( geneIdStr != lastGeneIdStr ) {
      ++geneId;
    }
    string strnd = get_tab_elem(buff,3);
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
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue: ["<<filename<<"]"<<std::endl;
    }
    lastGeneIdStr = geneIdStr;
  }
  if(!debug_cnt) {
    throw std::invalid_argument(filename + " : unable to read file");
  }
  return newGenePredDataSet;
}
