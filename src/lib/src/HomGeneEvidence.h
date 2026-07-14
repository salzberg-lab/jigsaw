//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef HOMGENEEVIDENCE_H
#define HOMGENEEVIDENCE_H

#include "EvidenceType.h"
#include "DataSetAbstractProduct.h"
#include <iomanip>
#include <math.h>

/**
 * Defines the behavior of homology evidence
 * in the feature space
 */
class HomGeneEvidence : public EvidenceType {
 public:
  HomGeneEvidence(const string&, int, double); 
  virtual ~HomGeneEvidence();

  virtual void finish();
  unsigned getNum() const;
  virtual void update(const Prediction* pred, int end5, int end3, dsu::Strand_t); 
  void setVals(const EvidenceType& et); 
  double eval() const; 
  virtual EvidenceType* dup() const; 
  bool train_print(std::ostream& os) const; 
  bool read(std::istream& is); 
};

class HomGeneEvidenceBndry : public HomGeneEvidence {
 public:
  HomGeneEvidenceBndry(const string&, int, double, const string&); 
  static bool isEnd5Splice(const DnaStr&, int, dsu::Strand_t); 
  static bool isEnd3Splice(const DnaStr&, int, dsu::Strand_t); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 

 protected:
  const string& _theSeq;
};

class HomGeneEvidenceBndryAcc : public HomGeneEvidenceBndry {
 public:
  HomGeneEvidenceBndryAcc(const string&, int, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomGeneEvidenceBndryDon : public HomGeneEvidenceBndry {
 public:
  HomGeneEvidenceBndryDon(const string&, int, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomGeneEvidenceBndryStart : public HomGeneEvidenceBndry {
 public:
  HomGeneEvidenceBndryStart(const string&, int, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomGeneEvidenceBndryStop : public HomGeneEvidenceBndry {
 public:
  HomGeneEvidenceBndryStop(const string&, int, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};

class HomGeneEvidenceIntron : public HomGeneEvidenceBndry {
 public:
  HomGeneEvidenceIntron(const string&, int, double, const string&); 
  virtual void update(const Prediction*, int, int, dsu::Strand_t); 
  virtual EvidenceType* dup() const; 
};


#endif
