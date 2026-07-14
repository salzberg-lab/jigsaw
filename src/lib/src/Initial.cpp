//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Initial.h"
#include <iostream>

Initial::Initial(){}
Initial Initial::_initial;

void Initial::print(std::ostream& os, const string& str) const {
  os<<str<<getTypeStr()<<": ";
}


static const string str = "Initial";
const string& Initial::getTypeStr() const { return str; }
