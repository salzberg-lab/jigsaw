//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEPREDEVIDENCE_H
#define GENEPREDEVIDENCE_H

#include "EvidenceType.h"
#include "HomologyEvidence.h" // is this the way?
#include "DataSetAbstractProduct.h"
#include "GenePrediction.h"
#include <ctype.h> //tolower

/**
 * Defines the behavior of gene prediction evidence as a feature
 */
class GenePredEvidence : public EvidenceType {
 public:
  GenePredEvidence(const string&, int, double);
  ~GenePredEvidence();

  virtual EvidenceType* dup() const; 
  virtual void update(const Prediction*, int, int, dsu::Strand_t);

  unsigned getNum() const { return 1; }
  //double eval() const; 
  bool train_print(std::ostream& os) const;
  bool read(std::istream&); 

};

class GenePredEvidenceBndry : public GenePredEvidence {
 public:
  GenePredEvidenceBndry(const string&, int, double);
  virtual EvidenceType* dup() const; 
  virtual void update(const Prediction*, int, int, dsu::Strand_t );
};

class GenePredEvidenceBndryAcc : public GenePredEvidenceBndry {
 public:
  GenePredEvidenceBndryAcc(const string& str, int idx, double val) :
    GenePredEvidenceBndry(str,idx,val) { }
  virtual EvidenceType* dup() const;
  virtual void update(const Prediction*, int, int , dsu::Strand_t);
};

class GenePredEvidenceBndryDon : public GenePredEvidenceBndry {
 public:
  GenePredEvidenceBndryDon(const string& str, int idx, double val) :
    GenePredEvidenceBndry(str,idx,val) { }
  virtual EvidenceType* dup() const;
  virtual void update(const Prediction*, int, int, dsu::Strand_t );
};

class GenePredEvidenceBndryStart : public GenePredEvidenceBndry {
 public:
  GenePredEvidenceBndryStart(const string& str, int idx, double val) :
    GenePredEvidenceBndry(str,idx,val) { }
  virtual EvidenceType* dup() const;
  virtual void update(const Prediction*, int, int, dsu::Strand_t );
};

class GenePredEvidenceOverlapStart : public GenePredEvidenceBndry {
 public:
  GenePredEvidenceOverlapStart(const string& str, int idx, double val) :
    GenePredEvidenceBndry(str,idx,val) { }
  virtual EvidenceType* dup() const;
  virtual void update(const Prediction*, int, int, dsu::Strand_t );
};

class GenePredEvidenceBndryStop : public GenePredEvidenceBndry {
 public:
  GenePredEvidenceBndryStop(const string& str, int idx, double val) :
    GenePredEvidenceBndry(str,idx,val) { }
  virtual EvidenceType* dup() const;
  virtual void update(const Prediction*, int, int, dsu::Strand_t );
};

class GenePredEvidenceIntron : public GenePredEvidenceBndry {
 public:
  GenePredEvidenceIntron(const string& str, int idx, double val) :
    GenePredEvidenceBndry(str,idx,val) { }
  virtual EvidenceType* dup() const;
  virtual void update(const Prediction*, int, int, dsu::Strand_t );
};

#endif
