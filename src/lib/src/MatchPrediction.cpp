//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "MatchPrediction.h"


void MatchPrediction::print(std::ostream& os, const std::string& str) const {
  GenePrediction::print(os,str);
  //os<<" "<<getPcntId();
  os<<" "<<getPcntSim();
}
