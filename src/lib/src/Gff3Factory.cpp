//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Gff3Factory.h"
#include "GenePredDataSet.h"
#include "SpliceSiteDataSet.h"
#include "HomologyDataSet.h"
#include "MatchPrediction.h"
#include "genemodels.h"
#include "DnaStr.h"
#include <fstream>
#include <stdexcept>
#include <list>
#include <map>

using std::cout;
using std::endl;

#if 0
static string 
conv2exon(const string& str) {
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
#endif

static string
getParentId(const string& comment) {
#if 0
  const string::size_type start = comment.find("Name=")+5;
  assert(start < comment.length());
  string::size_type stop = comment.find(";",start);
  if( string::npos == stop ) {
    stop = comment.length()-1;
  }
  assert(stop < comment.length());
  const string gid = comment.substr(start,stop-start);
  if( string::npos != gid.find(",") ) {
    throw std::runtime_error("Multiple parents are currently not handled, please assign one parent per feature");
  }
#endif
  const string::size_type start = comment.find("Parent=")+7;
  if( start == string::npos ) {
    throw std::runtime_error(string("can not find parent id on ") + comment);
  }
  assert(start < comment.length());
  string::size_type stop = comment.find(";",start);
  if( string::npos == stop ) {
    stop = comment.length()-1;
  }
  assert(stop < comment.length());
  const string gid = comment.substr(start,stop-start+1);
  if( string::npos != gid.find(",") ) {
    throw std::runtime_error("Multiple parents are currently not handled, please assign one parent per feature");
  }
  return gid;
}

#if 0
static string
getExonId(const string& comment) {
  const string::size_type start = comment.find("Note=")+5;
  assert(start < comment.length());
  string::size_type stop = comment.find(";",start);
  if( string::npos == stop ) {
    stop = comment.length()-1;
  }
  assert(stop < comment.length());
  const string gid = comment.substr(start,stop-start);
  return conv2exon(gid);
}
#endif

static string get_tab_elem(const string& str, unsigned ith_tab) {
  string::size_type pos=0,stop=string::npos; 
  if( pos >= str.size() ) {
    std::cerr<<"Gff3Factory:get_tab_elem failed"<<std::endl;
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

typedef std::map<string, int> local_hmap;

DataSetAbstractProduct* 
Gff3Factory::createDataSetAbstractProduct(const string& filename, const string& file_id, const string& evidence_type, const string& seq) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");
  if(strcmp(evidence_type.c_str(),"spliceprediction")==0) {
    throw std::invalid_argument("Gff3Factory error, spliceprediction format not available");
  } else if( (strcmp(evidence_type.c_str(),"homology")==0) || (strcmp(evidence_type.c_str(),"homgene")==0) ) {
    return gffHomologyRead(ifs,file_id,seq);
  }
  local_hmap htable;
  const int buff_size = 20096;
  char buff[buff_size];
  SeqLoc sl(dsu::eEither,-1,-1);
  // create a new list of genes
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,file_id);
  int debug_cnt = 0, geneId = 0, gcnt = 0;
  string gffid = file_id;
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '#' || *buff == '\n' || !*buff ) {
      continue;
    }
    const string exon_type = get_tab_elem(buff,2);
    if( !dsu::ignoreCaseEqual("CDS",exon_type) ) { //&& !dsu::ignoreCaseEqual("exon",exon_type) ) {
      continue;
    }
    const int end5 = atoi(get_tab_elem(buff,3).c_str());
    const int end3 = atoi(get_tab_elem(buff,4).c_str());
    const bool validEnd5 = (end5 >= 1 && end5 < (signed)seq.size());
    const bool validEnd3 = (end3 >= 1 && end3 < (signed)seq.size());
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
    const string buff_five = get_tab_elem(buff,5);
    double score = 1.0;
    if( buff_five != string(".") ) {
      score = atof(buff_five.c_str());
    }
    const string geneIdStr = get_tab_elem(buff,8);
    const string gstr = getParentId(geneIdStr);// + string("-padwiththis");
    const string exon_type_id = "internal";
    const local_hmap::const_iterator match = htable.find(gstr);
    if( match == htable.end() ) {
      ++gcnt;
      geneId = htable[gstr] = gcnt;
    } else {
      geneId = htable[gstr];
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
      newGenePredDataSet->insert(gffid,geneId,theSl,exon_type_id,score);
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue: ["<<filename<<"]"<<std::endl;
    }
    ++debug_cnt;
  }
  if(!debug_cnt) {
    throw std::invalid_argument(filename + " : unable to read file");
  }
  newGenePredDataSet->label_exons(seq);
  return newGenePredDataSet;
}


DataSetAbstractProduct* Gff3Factory::gffHomologyRead(std::ifstream& ifs, const string& type, const string& seq) {
  int debug_cnt = 0;
  HomologyDataSet::TranscriptHash* thPtr = new HomologyDataSet::TranscriptHash;
  HomologyDataSet* newHomologyDataSet = new HomologyDataSet(SeqLoc(dsu::eEither,-1,-1),type,thPtr);
  const int buff_size = 2024;
  char buff[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    if( *buff == '#' || *buff == '\n' ) 
      continue;
    ++debug_cnt;
    string exon_type = get_tab_elem(buff,2);
    if( !dsu::ignoreCaseEqual("CDS",exon_type) && !dsu::ignoreCaseEqual("exon",exon_type) ) {
      continue;
    }
    const string geneIdStr = get_tab_elem(buff,8);
    const string accession = getParentId(geneIdStr) + string("-") + exon_type;
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
  std::list<HomologyDataSet::iterator> del_lst;
  HomologyDataSet::iterator iter = newHomologyDataSet->begin();
  const HomologyDataSet::iterator stop = newHomologyDataSet->end();
  for(; iter != stop; ++iter) {
    HomologyDataSet* al_chain1 = dynamic_cast<HomologyDataSet*>(*iter);
    if( dsu::ignoreCaseEqual("exon",al_chain1->getDataSource()) ) {
      HomologyDataSet::iterator iter2 = newHomologyDataSet->begin();
      for(; iter2 != stop; ++iter2 ) {
	if( iter != iter2 ) {
	  const HomologyDataSet* al_chain2 = dynamic_cast<HomologyDataSet*>(*iter2);
	  if( isOverlap(*al_chain1,*al_chain2) && dsu::ignoreCaseEqual("CDS",al_chain2->getDataSource()) ) {
	    del_lst.push_back(iter);
	    break;
	  }
	}
      }
    }
  }
  std::list<HomologyDataSet::iterator>::iterator diter = del_lst.begin();
  std::list<HomologyDataSet::iterator>::iterator dstop = del_lst.end();
  for(; diter != dstop; ++diter) {
    Prediction* val = **diter;
    val->del();
    delete val;
    newHomologyDataSet->erase(*diter);
    
  }
  return newHomologyDataSet;
}
