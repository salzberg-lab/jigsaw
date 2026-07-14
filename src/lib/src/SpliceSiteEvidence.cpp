//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "SpliceSiteEvidence.h"
#include "DataSetAbstractProduct.h"
#include "SpliceSite.h"
#include <math.h>


/**
 * defines the behavior of splice site evidence as a feature
 *
 */
  SpliceSiteEvidence::SpliceSiteEvidence(const string& str, int idx, double wght) : 
    EvidenceType(str,idx,0,wght) { }

  SpliceSiteEvidence::~SpliceSiteEvidence() { }

  static inline double squash(double val) {
    const double base = 1.2;
    return 1-(1.0/pow(base,val));
  }
  unsigned SpliceSiteEvidence::getNum() const { return 1; }
  void SpliceSiteEvidence::update(const Prediction* pred,int end5, int end3, dsu::Strand_t strnd) { 
    if( pred->getStrnd() != strnd ) {
      return;
    }
    const SpliceSite* ss = dynamic_cast<const SpliceSite*>(pred);
    assert(ss);
    if(end5 != -1) {
      if( ss->getStrnd() == dsu::eNeg && ss->isDonor() ) {
	getValRef() = squash(ss->getScore().getVal());
      } else if(ss->getStrnd() == dsu::ePos && ss->isAcceptor()) {
	getValRef() = squash(ss->getScore().getVal());
      }
    } else if(end3 != -1) {
      if(ss->getStrnd() == dsu::eNeg && ss->isAcceptor()) {
	getValRef() = squash(ss->getScore().getVal());
      } else if(ss->getStrnd() == dsu::ePos && ss->isDonor()) {
	getValRef() = squash(ss->getScore().getVal());
      }
    }
  }

  EvidenceType* SpliceSiteEvidence::dup() const {
    return new SpliceSiteEvidence(getStr(),getIdx(),getWght());
  }
  bool SpliceSiteEvidence::train_print(std::ostream& os) const {
    os<<getVal(); 
    return true;
  }
  bool SpliceSiteEvidence::read(std::istream& is) {
    return (is>>getValRef()) ? true: false;
  }

  SpliceSiteAcc::SpliceSiteAcc(const string& str, int idx, double wght) : 
    SpliceSiteEvidence(str,idx,wght) { }

  EvidenceType* SpliceSiteAcc::dup() const {
    return new SpliceSiteAcc(getStr(),getIdx(),getWght());
  }

  void SpliceSiteAcc::update(const Prediction* pred,int end5, int end3, dsu::Strand_t strnd) { 
    if( pred->getStrnd() != strnd ) {
      return;
    }
    const SpliceSite* ss = dynamic_cast<const SpliceSite*>(pred);
    assert(ss);
    if( ss->isAcceptor() )
      getValRef() = squash(ss->getScore().getVal());
  }

  SpliceSiteDon::SpliceSiteDon(const string& str, int idx, double wght) : 
    SpliceSiteEvidence(str,idx,wght) { }

  EvidenceType* SpliceSiteDon::dup() const {
    return new SpliceSiteDon(getStr(),getIdx(),getWght());
  }
void SpliceSiteDon::update(const Prediction* pred,int end5, int end3, dsu::Strand_t strnd) { 
    if( pred->getStrnd() != strnd ) {
      return;
    }
    const SpliceSite* ss = dynamic_cast<const SpliceSite*>(pred);
    assert(ss);
    if( ss->isDonor() )
      getValRef() = squash(ss->getScore().getVal());
  }
