//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Prediction.h"

void Prediction::print(std::ostream& os, const string& str) const {
   os<<str;
   if(getModelId() != -1)
      os<<"model id: "<<getModelId()<<" ";
   os<<getTypeStr()<<" "<<getDataSource()<<" ";
   SeqLoc::print(os);
   os<<" "<<_score;
}

AnnotationScore Prediction::getScore() const { 
  return _score;
}
