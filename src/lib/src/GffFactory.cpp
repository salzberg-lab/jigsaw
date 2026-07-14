//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "GffFactory.h"
#include "GenePredDataSet.h"
#include "SpliceSiteDataSet.h"
#include "HomologyDataSet.h"
#include "MatchPrediction.h"
#include "genemodels.h"
#include "DnaStr.h"
#include <fstream>
#include <stdexcept>

using std::cout;
using std::endl;

static DataSetAbstractProduct* gffSpliceRead(std::ifstream&,const string&, const string&);

static string conv2exon(const string& str) {
  if(strstr(str.c_str(),"final")) return "terminal"; 
  else if(strstr(str.c_str(),"term")) return "terminal"; 
  else if(strstr(str.c_str(),"initial")) return "initial"; 
  else if(strstr(str.c_str(),"internal")) return "internal"; 
  else if(strstr(str.c_str(),"single")) return "single"; 
  else if(strstr(str.c_str(),"first")) return "initial"; 
  else if(strstr(str.c_str(),"last")) return "initial"; 
  else if(strstr(str.c_str(),"start")) return "initial"; 
  else if(strstr(str.c_str(),"stop")) return "terminal"; 
  else return "internal";
}

static string get_tab_elem(const string& str, unsigned ith_tab) {
  string::size_type pos=0,stop=string::npos; 
  if( pos >= str.size() ) {
    std::cerr<<"GffFactory:get_tab_elem failed"<<std::endl;
    std::cerr<<"["<<str<<"]"<<std::endl;
    std::cerr<<pos<<" "<<ith_tab<<std::endl;
    return "";
  }
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

DataSetAbstractProduct* GffFactory::createDataSetAbstractProduct(const string& filename, const string& file_id, const string& evidence_type, const string& seq) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");
  if(strcmp(evidence_type.c_str(),"spliceprediction")==0) {
    return gffSpliceRead(ifs,file_id,seq);
  } else if( (strcmp(evidence_type.c_str(),"homology")==0) || (strcmp(evidence_type.c_str(),"homgene")==0) ) {
    return gffHomologyRead(ifs,file_id,seq);
  }
  const int buff_size = 4096;
  char buff[buff_size];
  SeqLoc sl(dsu::eEither,-1,-1);
  // create a new list of genes
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,file_id);
  int debug_cnt = 0, geneId = 0;
  string geneIdStr,lastGeneIdStr;
  string gffid = file_id;
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '#' || *buff == '\n' || !*buff ) 
      continue;
    ++debug_cnt;
    string exon_type = get_tab_elem(buff,2);
    exon_type = conv2exon(exon_type);
    int end5 = atoi(get_tab_elem(buff,3).c_str());
    int end3 = atoi(get_tab_elem(buff,4).c_str());
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
    string buff_five = get_tab_elem(buff,5);
    double score = 1.0;
    if( buff_five != string(".") ) {
      score = atof(buff_five.c_str());
    }
    string geneIdStr = get_tab_elem(buff,8);
    if( geneIdStr != lastGeneIdStr ) {
      ++geneId;
    }
    string strnd = get_tab_elem(buff,6);
    dsu::Strand_t parseStrnd = dsu::eEither; 
    if( strnd == "-" ) {
      parseStrnd = dsu::eNeg;
    } else {
      parseStrnd = dsu::ePos;
    }
    try {
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


DataSetAbstractProduct* GffFactory::gffHomologyRead(std::ifstream& ifs, const string& type, const string& seq) {
  int debug_cnt = 0;
  HomologyDataSet::TranscriptHash* thPtr = new HomologyDataSet::TranscriptHash;
  HomologyDataSet* newHomologyDataSet = new HomologyDataSet(SeqLoc(dsu::eEither,-1,-1),type,thPtr);
  const int buff_size = 2024;
  char buff[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '#' || *buff == '\n' ) 
      continue;
    ++debug_cnt;
    //string accession = get_tab_elem(buff,1)+ string(" ") + get_tab_elem(buff,1);
    const string accession = get_tab_elem(buff,8);
    const int end5 = atoi(get_tab_elem(buff,3).c_str());
    const int end3 = atoi(get_tab_elem(buff,4).c_str());
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
    double pcid = 100.0;
    const string pscore = get_tab_elem(buff,5);
    if( pscore != string(".") ) {
      pcid = atof(get_tab_elem(buff,5).c_str());
    }
    const string strnd = get_tab_elem(buff,6);
    dsu::Strand_t parseStrnd = dsu::eEither; 
    if( strnd == "-" ) {
      parseStrnd = dsu::eNeg;
    } else {
      parseStrnd = dsu::ePos;
    }
    try {
      assert(end5 >= 0 && end3 >= 0);
      SeqLoc sl(parseStrnd,end5,end3);
      SeqLoc rsl( dsu::eEither, -1, -1 );
      const double pcsim = pcid, lscore = pcid, native_score = pcid;
      MatchPrediction* mp = new MatchPrediction(sl,native_score,-1,accession,pcid,pcsim,lscore,rsl,Internal::getInstance());
      newHomologyDataSet->insert(mp,type);
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue..."<<std::endl;
    }
  }
  if(!debug_cnt) 
    throw std::logic_error(type + " unable to get any data from this file");
  return newHomologyDataSet;
}


DataSetAbstractProduct* gffSpliceRead(std::ifstream& ifs, const string& type, const string& seq) {
  SeqLoc sl(dsu::eEither,-1,-1);
  SpliceSiteDataSet* newSpliceSiteDataSet = new SpliceSiteDataSet(sl,-1,type);
  int debug_cnt = 0;
  const int buff_size = 4096;
  char buff[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '#' || *buff == '\n' ) 
      continue;
    string gffid = get_tab_elem(buff,0)+ string(" ") + get_tab_elem(buff,1);
    string splice_type = get_tab_elem(buff,2);
    int end5 = atoi(get_tab_elem(buff,3).c_str());
    int end3 = atoi(get_tab_elem(buff,4).c_str());
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
    double score = atof(get_tab_elem(buff,5).c_str());
    string strnd = get_tab_elem(buff,6);
    if(strnd == "-") {
      std::swap(end5,end3);
    }
    ++debug_cnt;
    try {
      assert(end5 >= 0 && end3 >= 0);
      newSpliceSiteDataSet->insert(type, end5,end3,splice_type,score);
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue..."<<std::endl;
    }
  }
  if(!debug_cnt) 
    throw std::logic_error(type + " unable to get any data from this file");
    
  //newSpliceSiteDataSet->begin();
  return newSpliceSiteDataSet;
}
