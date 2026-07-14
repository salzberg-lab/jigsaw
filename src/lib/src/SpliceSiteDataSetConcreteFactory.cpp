//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "SpliceSiteDataSetConcreteFactory.h"
#include "SpliceSiteDataSet.h"
#include "DnaStr.h"
#include <fstream>
#include <stdexcept>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;

static int cntSpaces(const char* str) {
  int cnt = 0;
  for(const char* strPtr = str; *strPtr != '\0'; ++strPtr) {
    if(*strPtr == ' ')
      ++cnt;
  }
  return cnt;
}

DataSetAbstractProduct* 
SpliceSiteDataSetConcreteFactory::createDataSetAbstractProduct(const string& filename, const string& type, const string& /* evidence type */, const string& seq) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");
  
  string strnd, splice_type;
  int end5, end3;
  double score = 1.0;
  SeqLoc sl(dsu::eEither,-1,-1);
  SpliceSiteDataSet* newSpliceSiteDataSet = new SpliceSiteDataSet(sl,-1,type);
  int debug_cnt = 0;
  const int buff_size = 2024;
  char buff[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    int num_spaces = cntSpaces(buff);
    std::istringstream istrm(buff);
    if(num_spaces == 2 ) {
      istrm>>end5>>end3>>splice_type;
    } else if(num_spaces == 4 ) {
      string dumb;
      istrm>>end5>>end3>>score>>dumb>>splice_type;
    } else {
      istrm>>end5>>end3>>splice_type>>score;
    }
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
    try {
      assert(end5 >= 0 && end3 >= 0);
      newSpliceSiteDataSet->insert(type, end5,end3,splice_type,score);
    } catch(std::exception& e) {
      cerr<<e.what()<<endl;
      cerr<<"attempting to continue..."<<endl;
    }
  }
  if(!debug_cnt) 
    throw std::logic_error(filename + " unable to get any data from this file");
    
  newSpliceSiteDataSet->begin();
  return newSpliceSiteDataSet;
}
