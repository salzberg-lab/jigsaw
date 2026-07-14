//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "BtabFactory.h"
#include "HomologyDataSet.h"
#include "MatchPrediction.h"
#include "Internal.h"
#include <fstream>
#include <stdexcept>

using std::cout;
using std::endl;


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
/**
 **
#  0 : 51346 fixath1     #### asmbl_id querydatabase
#  1 : Jan 16 2000       #### date
#  2 : 99758             #### query sequence length
#  3 : gap2              #### method
#  4 : /usr/local/db/GUDB.arab    #### searchdatabase
#  5 : arab|TC56431      #### accession
#  6 : 10824             #### end5
#  7 : 10718             #### end3
#  8 : 529               #### match left
#  9 : 635               #### match right
#  10 : 68.224299        #### per_id
#  11 : 68.224299        #### per_sim
#  12 : 78               #### score
#  13 : -1               #### 
#  14 : -1               #### 
#  15 :                  #### description
#  16 : -1               #### 
#  17 : Minus            #### strand 
#  18 : 1531             #### match length
**/

DataSetAbstractProduct* 
BtabFactory::createDataSetAbstractProduct(const string& filename, const string& type, const string& /* evidence type*/, const string& seq) {
  std::ifstream ifs(filename.c_str());
  if(!ifs)
    throw std::invalid_argument(filename + " is not a valid input file ");
  
  int debug_cnt = 0;
  string strnd, accession;
  int relEnd5=-1, relEnd3=-1, native_score = 0;
  HomologyDataSet* newHomologyDataSet = new HomologyDataSet(SeqLoc(dsu::eEither,-1,-1),filename,NULL);
  const int buff_size = 2024;
  char buff[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    const string accession = get_tab_elem(buff,5);
    int end5 = atoi(get_tab_elem(buff,6).c_str());
    int end3 = atoi(get_tab_elem(buff,7).c_str());
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
    const double percid = atof(get_tab_elem(buff,10).c_str());
    const double percsim = atof(get_tab_elem(buff,11).c_str());
    const double lscore = atof(get_tab_elem(buff,12).c_str());
    ++debug_cnt;
    try {
      dsu::Strand_t theStrand = dsu::eEither;
      if(end5 > end3) {
	std::swap(end5,end3);
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
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue..."<<std::endl;
    }
  }
  if(!debug_cnt) 
    throw std::logic_error(filename + " unable to get any data from this file");

  return newHomologyDataSet;
}
