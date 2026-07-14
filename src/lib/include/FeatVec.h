//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef FEATVEC_H
#define FEATVEC_H

#include "EvidenceType.h"
#include <vector>

using std::vector;

/**
 * represents a vector of EvidenceType* 
 * and is used to store a heterogenous collection
 * evidence types.
 *
 */

class FeatVec : public vector<EvidenceType*> {
 public:
  /**
   * feature vectors should be used as a fixed size vector
   * @param int size - number of elements in the feature fector
   */
  FeatVec(int size) : vector<EvidenceType*>(size) { }
  FeatVec(const FeatVec&);

  /*
   * calls delete for on each element
   */
  ~FeatVec();

  /**
   * creates an orthonormal vector
   */
  void normalize();

  double voteEvalOpp() const;
  double voteEval() const;
  double sumEval() const;
  double prodEval() const;
  double avgEval() const;

  /**
   * return the number of feature dimensions for the feature vector
   */
  int getNumDim() const;

  /**
   * convert feature vector to a array of floats
   */
  float* toFloatArr(float*,unsigned) const;

  /**
   * some EvidenceType objects might have a normalization feature
   * that needs to be run
   */
  void elem_normalize();

  /**
   * set the values of the feature vector
   * @param - const FeatVec& input values to set 
   */
  void setVals(const FeatVec&);

  /**
   * print the values
   */
  std::ostream& val_print(std::ostream&) const;

  /**
   * return the sum of the difference between each matching element of two vectors
   * return (\Sum xi-yi)
   */
  friend double operator-(const FeatVec&, const FeatVec&);

  void finish();

};
#endif // FEATVEC_H
