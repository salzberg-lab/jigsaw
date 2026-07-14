//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef _ANSCHECK_H
#define _ANSCHECK_H

#include "DataSetAbstractProduct.h"
#include "RegionPrediction.h"
#include "dsu.h"
#include <vector>

class AnsCheck {
 public:
  static bool checkCuration(const RegionPrediction&,dsu::Strand_t);
  static void onBorder(const DataSetAbstractProduct&,const RegionPrediction&, dsu::Strand_t, vector<bool>&);
};

#endif //_ANSCHECK_H
