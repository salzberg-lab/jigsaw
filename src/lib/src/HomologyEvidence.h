//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef HOMOLOGYEVIDENCE_H
#define HOMOLOGYEVIDENCE_H

#include "EvidenceType.h"
#include "DataSetAbstractProduct.h"
#include <iomanip>
#include <math.h>

/**
 * Defines the behavior of homology evidence
 * in the feature space
 */
class HomologyEvidence : public EvidenceType {
 public:
  HomologyEvidence(const string&, int, double, double, double); 
  virtual ~HomologyEvidence();

  virtual void finish();
  unsigned getNum() const;
  virtual void updateSum(double val); 
  virtual void update(const Prediction* pred, int end5, int end3, dsu::Strand_t); 
  void setVals(const EvidenceType& et); 
  //double eval() const; 
  virtual EvidenceType* dup() const; 
  bool train_print(std::ostream& os) const; 
  bool read(std::istream& is); 
  double getVal(int idx) const; 
 protected:
  double _sumVal, _numWght, _norm;
};

class HomologyEvidenceBndry : public HomologyEvidence {
 public:
  HomologyEvidenceBndry(const string&, int, double, double, double, const string&); 
  static bool isEnd5Splice(const DnaStr&, int, dsu::Strand_t); 
  static bool isEnd3Splice(const DnaStr&, int, dsu::Strand_t); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 

 protected:
  const string& _theSeq;
};

class HomologyEvidenceBndryAcc : public HomologyEvidenceBndry {
 public:
  HomologyEvidenceBndryAcc(const string&, int, double, double, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomologyEvidenceBndryDon : public HomologyEvidenceBndry {
 public:
  HomologyEvidenceBndryDon(const string&, int, double, double, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomologyEvidenceBndryStart : public HomologyEvidenceBndry {
 public:
  HomologyEvidenceBndryStart(const string&, int, double, double, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomologyEvidenceBndryStop : public HomologyEvidenceBndry {
 public:
  HomologyEvidenceBndryStop(const string&, int, double, double, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomologyEvidenceIntron : public HomologyEvidenceBndry {
 public:
  HomologyEvidenceIntron(const string&, int, double, double, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
  void finish();
 private:
  int _intronPred, _nonIntronPred;
};


#endif
