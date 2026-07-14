//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "AnnotationScore.h"
#include <iostream>
#include <iomanip>

std::ostream& operator<<(std::ostream& os, const AnnotationScore& score) {
  os<<std::setprecision(10)<<score.getVal();
  return os;
}

std::istream& operator>>(std::istream& is, AnnotationScore& score) {
  is>>score._val;
  return is;
}
