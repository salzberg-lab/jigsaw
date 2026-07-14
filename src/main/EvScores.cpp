//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "EvScores.h"
#include "AnnotationItem.h"
#include "genemodels.h"
#include <iostream>
#include <vector>

EvScores::EvScores() {
  for(unsigned i = 0; i < eNumBnd; ++i) {
    (*this)[i] = AnnotationScore::badVal();
  }
}

EvScores::EvScores(const EvScores& evs) {
  for(unsigned i = 0; i < eNumBnd; ++i) {
    (*this)[i] = evs[i];
  }
}

EvScores& EvScores::operator=(const EvScores& evs) {
  for(unsigned i = 0; i < eNumBnd; ++i) {
    (*this)[i] = evs[i];
  }
  return *this;
}

std::ostream& operator<<(std::ostream& os, const EvScores& es) {
  for(unsigned i = 0; i < EvScores::eNumBnd; ++i) {
    os<<(es[i])<<" ";
  }
  return os;
}
