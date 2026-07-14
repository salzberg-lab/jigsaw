//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "GenePredDataSetConcreteFactory.h"
#include "GenePredDataSet.h"
#include "genemodels.h"
#include "DnaStr.h"
#include <fstream>
#include <stdexcept>
#include <sstream>

using std::list;
using std::endl;
using std::cout;

static DataSetAbstractProduct* createDataSetAbstractProduct(const string& , const string& );
static DataSetAbstractProduct* createDataSetAbstractProduct2(const string& , const string& );

static bool useCurationFormat(char* buff) {
  int cnt = 0;
  for(char* ptr = buff; *ptr && *ptr != '\n'; ++ptr) {
    if(*ptr == ' ') ++cnt;
  }
  if(cnt==2) return true; 
  return false;
}

static bool noExonType(const string& buff) {
  bool status = true;
  string::size_type first = buff.find(' ')+1;
  string::size_type secnd = buff.find(' ',first);
  if( string::npos != first && string::npos != secnd ) {
    string sstr = buff.substr(first,secnd-first);
    if(
       sstr == "Internal" ||
       sstr == "internal" ||
       sstr == "Terminal" ||
       sstr == "terminal" ||
       sstr == "Initial" ||
       sstr == "initial" ||
       sstr == "Single" ||
       sstr == "single" ) {
      status =false;
    } else {
      status = true;
    }
  }
  return status;
}


DataSetAbstractProduct* GenePredDataSetConcreteFactory::createDataSetAbstractProduct(const string& filename, const string& type, const string& evtype, const string& seq){
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");

  char buff[500];
  ifs.getline(buff,500);
  if( useCurationFormat( buff) ) {
    return ::createDataSetAbstractProduct(filename,type);
  }
  if( noExonType(buff) ) {
    return ::createDataSetAbstractProduct2(filename,type);
  }
  ifs.seekg(0);
  string strnd;
  int geneId, end5, end3;
  SeqLoc sl(dsu::eEither,-1,-1);
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,type);
  int debug_cnt = 0;
  double prob;
  string dumbtype;
  while( ifs>>geneId>>dumbtype>>end5>>end3>>prob ) {
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
    ++debug_cnt;
    dsu::Strand_t parseStrnd = dsu::eEither; 
    if( end5 > end3 ) {
      std::swap(end5,end3);
      parseStrnd = dsu::eNeg;
    } else {
      parseStrnd = dsu::ePos;
    }
    try {
      assert(end5 >= 0 && end3 >= 0);
      SeqLoc theSl(parseStrnd,end5,end3);
      newGenePredDataSet->insert(type,geneId,theSl,dumbtype,prob);
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue..."<<std::endl;
    }
  }
  if(!debug_cnt) {
    throw std::invalid_argument(filename + " : unable to read file");
  }
  return newGenePredDataSet;
}


static const string& 
getType(const list<int>::const_iterator& curr, const list<int>::const_iterator& last, 
	const list<int>::const_iterator& end, dsu::Strand_t strnd) {
  int last_id = last==end?-1:*last;
  int middle = *curr;
  list<int>::const_iterator tmp = curr;
  ++tmp;
  int third = tmp == end ? -1:*tmp;
  if( last_id == middle ) {
     if(third == -1 || third != middle) 
       return strnd==dsu::ePos ? Terminal::getInstance().getTypeStr() : Initial::getInstance().getTypeStr();
      else
	return Internal::getInstance().getTypeStr();
  } else if(last_id == -1) {
    if(third == -1 || third != middle) 
      return Single::getInstance().getTypeStr();
    else
      return strnd==dsu::ePos ? Initial::getInstance().getTypeStr() : Terminal::getInstance().getTypeStr();
  } else {
    if(third == -1 || third != middle)
      return Single::getInstance().getTypeStr();
    else
      return strnd==dsu::ePos ? Initial::getInstance().getTypeStr() : Terminal::getInstance().getTypeStr();
  }
}

static DataSetAbstractProduct* createDataSetAbstractProduct(const string& filename, const string& type) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");
  
  SeqLoc sl(dsu::eEither,-1,-1);
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,type);
  string label = type;

  list<int> geneIdList, end5List, end3List;
  int geneId, end5, end3, debug_cnt = 0;
  while( ifs>>geneId>>end5>>end3) {
    ++debug_cnt;
    geneIdList.push_back(geneId);
    end5List.push_back(end5);
    end3List.push_back(end3);
  }
  if(!debug_cnt) 
    throw std::logic_error(filename + " unable to get any data from this file");
  
  list<int>::const_iterator geneIdIterator = geneIdList.begin(), end5Iterator = end5List.begin();
  list<int>::const_iterator end3Iterator = end3List.begin(), last_iter = geneIdList.end();
  while( geneIdIterator != geneIdList.end() ) {
    // note there is a slight non-fatal bug here: even though the initial and terminal
    // exons are incorrectly assigned on the negative strand - they just need to be labeled
    // "not internal" to be processed correctly - may clean up at some point.
    end5 = *end5Iterator;
    end3 = *end3Iterator;
    geneId = *geneIdIterator; 
    dsu::Strand_t strnd = dsu::coord2strnd(end5,end3);
    if(end5 > end3) {
      int tmp = end5;
      end5 = end3;
      end3 = tmp;
    }
    const string& ptype = getType(geneIdIterator,last_iter,geneIdList.end(),strnd);
    try {
      SeqLoc theSl(strnd,end5,end3);
      newGenePredDataSet->insert(label,geneId,theSl,ptype,1.0);
      last_iter = geneIdIterator;
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue..."<<std::endl;
    }
    ++geneIdIterator;
    ++end5Iterator;
    ++end3Iterator;
  }
  return newGenePredDataSet;
}

static DataSetAbstractProduct* createDataSetAbstractProduct2(const string& filename, const string& type) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");
  
  SeqLoc sl(dsu::eEither,-1,-1);
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,type);
  string label = type;

  list<double> probLst;
  list<int> geneIdList, end5List, end3List;
  int geneId, end5, end3, debug_cnt = 0;
  //string strId, lastStrId;
  //const int buff_size = 2048;
  //char buff[buff_size];
  double prob;
  while( ifs>>geneId>>end5>>end3>>prob) {
    ++debug_cnt;
    geneIdList.push_back(geneId);
    end5List.push_back(end5);
    end3List.push_back(end3);
    probLst.push_back(prob);
  }
  if(!debug_cnt) 
    throw std::logic_error(filename + " unable to get any data from this file");
  
  list<double>::const_iterator pi = probLst.begin();
  const list<double>::const_iterator ps = probLst.end();
  list<int>::const_iterator geneIdIterator = geneIdList.begin(), end5Iterator = end5List.begin();
  list<int>::const_iterator end3Iterator = end3List.begin(), last_iter = geneIdList.end();
  while( geneIdIterator != geneIdList.end() ) {
    // note there is a slight non-fatal bug here: even though the initial and terminal
    // exons are incorrectly assigned on the negative strand - they just need to be labeled
    // "not internal" to be processed correctly - may clean up at some point.
    assert(pi != ps);
    prob = *pi;
    end5 = *end5Iterator;
    end3 = *end3Iterator;
    geneId = *geneIdIterator; 
    dsu::Strand_t strnd = dsu::coord2strnd(end5,end3);
    if(end5 > end3) {
      int tmp = end5;
      end5 = end3;
      end3 = tmp;
    }
    const string& ptype = getType(geneIdIterator,last_iter,geneIdList.end(),strnd);
    try {
      SeqLoc theSl(strnd,end5,end3);
      newGenePredDataSet->insert(label,geneId,theSl,ptype,prob);
      last_iter = geneIdIterator;
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue..."<<std::endl;
    }
    ++geneIdIterator;
    ++end5Iterator;
    ++end3Iterator;
    ++pi;
  }
  return newGenePredDataSet;
}
