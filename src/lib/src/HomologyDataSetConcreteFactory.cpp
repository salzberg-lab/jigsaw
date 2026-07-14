//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "HomologyDataSetConcreteFactory.h"
#include "HomologyDataSet.h"
#include "MatchPrediction.h"
#include "Internal.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

using std::invalid_argument;
using std::ifstream;
using std::istringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::swap;

static int cntSpaces(const char* str) {
  int cnt = 0;
  for(const char* strPtr = str; *strPtr != '\0'; ++strPtr) {
    if(*strPtr == ' ')
      ++cnt;
  }
  return cnt;
}

DataSetAbstractProduct* HomologyDataSetConcreteFactory::createDataSetAbstractProduct(const string& filename, const string& type, const string& /* evidence type*/, const string& seq) {
  ifstream ifs(filename.c_str());
  if(!ifs)
    throw invalid_argument(filename + " is not a valid input file ");
  
  int debug_cnt = 0;
  string strnd, accession;
  int end5, end3, relEnd5, relEnd3, native_score = 0;
  double percid, percsim;
  HomologyDataSet* newHomologyDataSet = new HomologyDataSet(SeqLoc(dsu::eEither,-1,-1),filename,NULL);
  const int buff_size = 2024;
  char buff[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    int num_spaces = cntSpaces(buff);
    istringstream istrm(buff);
    int lscore = 1;
    if(num_spaces == 6 ) {
      istrm>>end5>>end3>>relEnd5>>relEnd3>>accession>>percid>>percsim;
    } else {
      istrm>>end5>>end3>>relEnd5>>relEnd3>>accession>>percid>>percsim>>lscore;
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
      dsu::Strand_t theStrand = dsu::eEither;
      if(end5 > end3) {
	swap(end5,end3);
	theStrand = dsu::eNeg;
      } else {
	theStrand = dsu::ePos;
      }
      assert(end5 >= 0 && end3 >= 0);
      SeqLoc sl(theStrand,end5,end3);
      dsu::Strand_t relStrand = dsu::ePos;
      int left = relEnd5, right = relEnd3;
      if( relEnd5 > relEnd3 ) {
	left = relEnd3;
	right = relEnd5;
	relStrand = dsu::eNeg;
      }
      SeqLoc rsl( relStrand, left, right );
      MatchPrediction* mp = new MatchPrediction(sl,native_score,-1,accession,percid,percsim,lscore,rsl,Internal::getInstance());
      newHomologyDataSet->insert(mp,type);
      //** TODO
      //** stop making the genepredata extra indirection
      //** but you have to correct the checkForMatchPrediction function
      //** and then retest
    } catch(std::exception& e) {
      cerr<<e.what()<<endl;
      cerr<<"attempting to continue..."<<endl;
    }
  }
  if(!debug_cnt) 
    throw std::logic_error(filename + " unable to get any data from this file");

  return newHomologyDataSet;
}
