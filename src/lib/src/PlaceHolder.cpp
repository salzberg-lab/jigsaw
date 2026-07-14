//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "PlaceHolder.h"

using std::endl;

const string PlaceHolder::_ph = "PlaceHolder";

void PlaceHolder::print(std::ostream& os, const std::string& str) const {
  Prediction::print(os,str);
  os<<std::endl<<"Points to this: "<<endl;
  _totalSequence->print(os,str);
}

