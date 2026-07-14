//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "AnnotationMatrix.h"
#include "CustomScore.h"
#include "AnnotationItem.h"
#include "SeqLoc.h"
#include "genemodels.h"
#include "ExonDistr.h"
#include "Options.h"
#include "Compare.h"
#include <stdexcept>

using std::cout;
using std::endl;

AnnotationMatrix::AnnotationMatrix() {
  const Options& opt = *Options::getInstance();
  _cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    _cmp = GtCompare::getInstance();

}

AnnotationMatrix::~AnnotationMatrix() {
  del();
}

void AnnotationMatrix::del() {
  AnnotationMatrix::const_iterator iter = end(), stop = end();
  while(iter != stop) {
    AnnotationItem* item = *iter;
    delete item;
    ++iter;
  }
  _data.clear();
}

void AnnotationMatrix::insert(AnnotationItem* item) {
  _data.push_back(item);
}

AnnotationMatrix::iterator AnnotationMatrix::begin() {
  return _data.begin();
}

AnnotationMatrix::iterator AnnotationMatrix::end() {
  return _data.end();
}

AnnotationMatrix::const_iterator AnnotationMatrix::begin() const {
  return _data.begin();
}

AnnotationMatrix::const_iterator AnnotationMatrix::end() const {
  return _data.end();
}
