//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Intron.h"
#include <iostream>

Intron Intron::_intron;

void Intron::print(std::ostream& os,const string& str) const {
  os<<str<<getTypeStr()<<": ";
}

static const string singleStr = "Intron";
const string& Intron::getTypeStr() const {
  return singleStr;
}
Intron::Intron(){}
