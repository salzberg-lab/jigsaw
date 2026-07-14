//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Terminal.h"
#include <iostream>

Terminal::Terminal(){}
Terminal Terminal::_terminal;

void Terminal::print(std::ostream& os,const string& str) const {
  os<<str<<getTypeStr()<<": ";
}

static const string terminalStr = "Terminal";
const string& Terminal::getTypeStr() const {
  return terminalStr;
}
