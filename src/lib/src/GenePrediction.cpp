//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "GenePrediction.h"
#include "GeneModelType.h"
#include "PlaceHolder.h"

using std::endl;

Prediction::DataStorageType* GenePrediction::makePlaceHolders() const {
  DataStorageType* retVal = new DataStorageType; 
  PlaceHolder* ph = new PlaceHolder(static_cast<const Prediction&>(*this), static_cast<const Prediction*>(this));
  retVal->push_back( ph );
  return retVal;
}

void GenePrediction::print(std::ostream& os, const string& str) const {
  Prediction::print(os," ");
}

const string& GenePrediction::getTypeStr() const {
  return _modelType.getTypeStr();
}

const int GenePrediction::getScoreIdx() const { 
  return _modelType.getScoreIdx();
}
