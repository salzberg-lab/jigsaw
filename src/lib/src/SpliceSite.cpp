//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "SpliceSite.h"
#include "genemodels.h"
#include "GenePrediction.h"
#include "PlaceHolder.h"
#include <iostream>

void SpliceSite::print(std::ostream& os, const string& str) const {
  os<<str<<_spliceType<<" : ";
  Prediction::print(os," ");
}

Prediction::DataStorageType* SpliceSite::makePlaceHolders() const {
  DataStorageType* retVal = new DataStorageType;
  Prediction* ph = new PlaceHolder(*this,this);
  retVal->push_back( ph );
  return retVal;
}

static const string str = "SpliceSite";
const string& SpliceSite::getTypeStr() const { return str; }

const string& SpliceSite::getSpliceType() const { return _spliceType; }

bool SpliceSite::isAcceptor() const { 
  return (getSpliceType()=="donor") ? false : true;
} 

bool SpliceSite::isDonor() const { 
  return (getSpliceType()=="donor") ? true : false;
} 
