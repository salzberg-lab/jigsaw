//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Igr.h"
#include <iostream>

Igr Igr::_igr;

void Igr::print(std::ostream& os,const string& str) const {
  os<<str<<getTypeStr()<<": ";
}

static const string singleStr = "Igr";
const string& Igr::getTypeStr() const {
  return singleStr;
}
Igr::Igr(){}
