//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "FeatVec.h"
#include <math.h>

FeatVec::FeatVec(const FeatVec& vec) : vector<EvidenceType*>(vec.size()) {
  FeatVec& ev_vec = *this;
  for(unsigned i = 0; i < size(); ++i) {
    ev_vec[i] = vec[i]->dup();
  }
}

FeatVec::~FeatVec() {
  FeatVec& ev_vec = *this;
  for(unsigned i = 0; i < ev_vec.size(); ++i)
    delete ev_vec[i];
}

double FeatVec::voteEval() const {
  double score = 0;
  const FeatVec& ev_vec = *this;
    // note this assumes confidence values (getVal()) range > 0 and weight range is [0,1]
  for(unsigned i = 0; i < ev_vec.size(); ++i) {
    //score += ev_vec[i]->getVal() > 0 ? ev_vec[i]->getWght() : 1-ev_vec[i]->getWght();;
    const double vote = ev_vec[i]->getVal();
    const double wght = ev_vec[i]->getWght();
    score += (vote*wght);
  }
  return score;
}

double FeatVec::voteEvalOpp() const {
  double score = 0;
  const FeatVec& ev_vec = *this;
  for(unsigned i = 0; i < ev_vec.size(); ++i) {
    // note this assumes confidence values (getVal()) range > 0 and weight range is [0,1]
    //score += ev_vec[i]->getVal() > 0 ? 1-ev_vec[i]->getWght() : ev_vec[i]->getWght();
    const double vote = 1-ev_vec[i]->getVal();
    const double wght = ev_vec[i]->getWght();
    score += (vote*wght);
  }
  return score;
}

double FeatVec::sumEval() const {
  double score = 0;
  const FeatVec& ev_vec = *this;
  for(unsigned i = 0; i < ev_vec.size(); ++i) 
    score += ev_vec[i]->eval();
  return score;
}

double FeatVec::avgEval() const {
  double score = 0;
  const FeatVec& ev_vec = *this;
  for(unsigned i = 0; i < ev_vec.size(); ++i) 
    score += ev_vec[i]->eval();
  return score/ev_vec.size();
}

double FeatVec::prodEval() const {
  double score = 1;
  const FeatVec& ev_vec = *this;
  for(unsigned i = 0; i < ev_vec.size(); ++i) 
    score *= ev_vec[i]->eval();
  return score;
}

int FeatVec::getNumDim() const {
  int num_d = 0;
  for(unsigned i = 0; i < size(); ++i) {
    const EvidenceType& ev = *((*this)[i]);
    num_d += ev.getNum();
  }
  return num_d;
}

float* FeatVec::toFloatArr(float* dimension, unsigned start) const {
  unsigned jmp = 0, mabs = start;
  for(unsigned i = 0; i < size(); ++i, mabs += jmp) {
    const EvidenceType& ev = *((*this)[i]);
    jmp = ev.getNum();
    for(unsigned j = 0; j < ev.getNum(); ++j) {
      dimension[mabs+j] = ev.getVal(j);
    }
  }
  return dimension;
}

void FeatVec::setVals(const FeatVec& vec) {
  for(unsigned i = 0; i < size(); ++i) {
    FeatVec& fv = *this;
    EvidenceType& et = *(fv[i]);
    et.setVals(*(vec[i]));
  }
}

void FeatVec::finish() {
  for(unsigned i = 0; i < size(); ++i) {
    this->operator[](i)->finish();
  }
}

void FeatVec::normalize() {
  double sum = 0;
  for(unsigned i = 0; i < size(); ++i) {
    const FeatVec& fv = *this;
    const EvidenceType& et = *(fv[i]);
    sum += et.getSumVal();
  }
  if(sum > 0) {
    for(unsigned i = 0; i < size(); ++i) {
      FeatVec& fv = *this;
      EvidenceType& et = *(fv[i]);
      et.normalize(sum);
    }
  }
}

std::ostream& FeatVec::val_print(std::ostream& os) const {
  const FeatVec& fv = *this;
  for(unsigned i = 0; i < fv.size(); ++i) {
    fv[i]->train_print(os);
    if(i < (size()-1))
      os<<" ";
  }
  return os;
}

double operator-(const FeatVec& afv, const FeatVec& bfv) {
  if(afv.size() != bfv.size())
     throw "finding distance between vectors of different demensions, not allowed";

  double dist = 0;
  for(unsigned i = 0; i < afv.size(); ++i) {
    double val = fabs(afv[i]->getVal()-bfv[i]->getVal());
    dist += val;
  }
  return dist;
}
