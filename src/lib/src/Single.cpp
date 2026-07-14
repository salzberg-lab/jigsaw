//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Single.h"
#include <iostream>

Single Single::_single;

void Single::print(std::ostream& os,const string& str) const {
  os<<str<<getTypeStr()<<": ";
  //os<<"warning for moment";
}

const string singleStr = "Single";
const string& Single::getTypeStr() const {
  return singleStr;
}
Single::Single(){}
