//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Internal.h"
#include <iostream>

Internal::Internal(){}
Internal Internal::_internal;

void Internal::print(std::ostream& os, const string& str) const {
  os<<str<<getTypeStr()<<": ";
}

const string internalStr = "Internal";
const string& Internal::getTypeStr() const {
  return internalStr;
}
