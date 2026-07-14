//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef _COMPARE_H
#define _COMPARE_H

#include "AnnotationScore.h"

/**
 ** class Compare
 ** 
 ** Allows the scoring function to be either a minimization or maximization.
 **
 ** Derived classes are LtCompare (less than or minimize) and 
 ** GtCompare (greater than or maximize)
 **/
class Compare {
 public:
  virtual bool compare(const AnnotationScore&,const AnnotationScore&) const = 0;
};

class LtCompare : public Compare {
 public:
  static const LtCompare* getInstance() {
    if(!_instance)
      _instance = new LtCompare();
    return _instance;
  }
  bool compare(const AnnotationScore& as1, const AnnotationScore& as2) const {
    if(as1.isValid() && as2.isValid())
       return as1<=as2;
    return as1.isValid();
  }
 private:
  static LtCompare* _instance;
};

class GtCompare : public Compare {
 public:
  static const GtCompare* getInstance() {
    if(!_instance)
      _instance = new GtCompare();
    return _instance;
  }
  bool compare(const AnnotationScore& as1, const AnnotationScore& as2) const {
    if(as1.isValid() && as2.isValid())
       return as1>=as2;
    else
       return as1.isValid();
  }
 private:
  static GtCompare* _instance;
};

#endif // _COMPARE_H
