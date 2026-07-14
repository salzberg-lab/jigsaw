//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef SPLICESITEEVIDENCE_H
#define SPLICESITEEVIDENCE_H

#include "EvidenceType.h"

/**
 * defines the behavior of splice site evidence as a feature
 *
 */
class SpliceSiteEvidence : public  EvidenceType {
 public:
  SpliceSiteEvidence(const string&, int, double);
  ~SpliceSiteEvidence();

  unsigned getNum() const;
  virtual void update(const Prediction* pred,int end5, int end3, dsu::Strand_t); 
  //double eval() const;
  EvidenceType* dup() const; 
  bool train_print(std::ostream& os) const; 
  bool read(std::istream& is); 

  //protected:
    //double _wght;
};

class SpliceSiteAcc : public  SpliceSiteEvidence {
 public:
  SpliceSiteAcc(const string& str, int idx, double wght); 

  EvidenceType* dup() const; 

  virtual void update(const Prediction* pred,int end5, int end3, dsu::Strand_t); 
};

class SpliceSiteDon : public  SpliceSiteEvidence {
 public:
  SpliceSiteDon(const string& str, int idx, double wght);

  EvidenceType* dup() const; 
  virtual void update(const Prediction* pred,int end5, int end3, dsu::Strand_t); 
};

#endif
