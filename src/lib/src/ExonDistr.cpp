//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "ExonDistr.h"
#include <assert.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <iostream>

using std::cerr;
using std::endl;
using std::sort;
using std::cout;

ExonDistr* ExonDistr::_instance = NULL;


static double calcMinFactor(double val) {
  assert(val >= 0 && val <= 1);
  if( val==0 || val==1 ) return val;
  double factor = 10;
  while( factor*val < 1 ) {
    factor *= 10;
  }
  return factor;
}

ExonDistr::ExonDistr(const string& filename) {
  if( filename.size() > 0 ) {
    std::ifstream ifs(filename.c_str());
    if( !ifs ) {
      cerr<<"warning, unable to open exon distr. file: "<<filename<<endl;
      cerr<<"ignoring...\n";
      return;
    }
    double val = -1;
    while(ifs>>val) {
      _pdistr.push_back(val);
    }
    std::vector<double> tmpvec = _pdistr;
    sort(tmpvec.begin(),tmpvec.end());
    double min = tmpvec[0];
    cout<<"min: "<<min<<endl;
    _minFactor = calcMinFactor(min);
    cout<<"minFactor: "<<_minFactor<<endl;
  }
}

double ExonDistr::getProb(int idx) const {
  // if 0, no values were provided
  // and we assume equa probability
  if( _pdistr.size() == 0 )
    return 0.0;

  if( idx < (int)_pdistr.size()) {
    double prob = 0.0;
    for(unsigned i = (unsigned)idx; i < _pdistr.size(); ++i) {
      prob += _pdistr[i];
    }
    assert(prob > 0);
    if( prob > 1 ) {
      cerr<<"shoot: "<<prob<<endl;
      assert(prob <= 1);
      exit(1);
    }
    return -log(_pdistr[idx]);
    //return -log(_pdistr[idx]);
  } else {
	 //The validity of this whole thing was nevery fully established
	  // henced the ugliness here...
    //unsigned li = _pdistr.size()-1;
    //return _pdistr[li];
    //double val = _minFactor > 0 ? ((double)1.0/(double)(_minFactor*idx)) : 0;
    //return val==0?100:-log(val);
    return 250;
  }
}

const ExonDistr& ExonDistr::getInstance() {
  return *_instance;
}

void ExonDistr::createInstance(const string& filename) {
  _instance = new ExonDistr(filename);
}
