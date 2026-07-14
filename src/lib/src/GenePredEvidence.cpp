//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "GenePredEvidence.h"
#include "genemodels.h"
#include "TypeManage.h"



GenePredEvidence::GenePredEvidence(const string& str, int idx, double wght) : 
  EvidenceType(str,idx,0,wght) { }

GenePredEvidence::~GenePredEvidence() { 
}

void GenePredEvidence::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  assert(end5 >= 0 && end3 >= 0);
  const GenePrediction* gp = TypeManage::checkForGenePrediction(pred,end5,end3);
  assert(gp);
  if( strnd == gp->getStrnd() && gp->getScoreIdx() != Intron::getInstance().getScoreIdx() ) {
    double val = gp->getScore().getVal();
    getValRef() = val;
  }
}

//double GenePredEvidence::eval() const {
  //return _wght*(getVal());
//}

EvidenceType* GenePredEvidence::dup() const {
  return new GenePredEvidence(getStr(),getIdx(),_wght);
}

bool GenePredEvidence::train_print(std::ostream& os) const {
  os<<getVal();
  return true;
}

bool GenePredEvidence::read(std::istream& is) {
  return is>>getValRef() ? true : false;
}

GenePredEvidenceBndry::GenePredEvidenceBndry(const string& str, int idx, double wght) :
  GenePredEvidence(str,idx,wght) {
}

EvidenceType* GenePredEvidenceBndry::dup() const {
    return new GenePredEvidenceBndry(getStr(),getIdx(),_wght);
}

void GenePredEvidenceBndry::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
    const GenePrediction* gp = TypeManage::checkForGenePrediction(pred,end5,end3);
    assert(gp);
    if(gp->getEnd5() == end5 && end5 != -1) getValRef() = gp->getScore().getVal();
    else if(gp->getEnd3() == end3 && end3 != -1) getValRef() = gp->getScore().getVal();
}

EvidenceType* GenePredEvidenceBndryAcc::dup() const {
  return new GenePredEvidenceBndryAcc(getStr(),getIdx(),_wght);
}

EvidenceType* GenePredEvidenceBndryDon::dup() const {
  return new GenePredEvidenceBndryDon(getStr(),getIdx(),_wght);
}

EvidenceType* GenePredEvidenceBndryStart::dup() const {
  return new GenePredEvidenceBndryStart(getStr(),getIdx(),_wght);
}

EvidenceType* GenePredEvidenceBndryStop::dup() const {
  return new GenePredEvidenceBndryStop(getStr(),getIdx(),_wght);
}

void GenePredEvidenceBndryAcc::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  assert(dp);
  const GenePrediction* gp = TypeManage::checkForGenePrediction(pred,end5,end3);
  assert(gp);
  if( gp->getStrnd() != strnd || gp->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
    
  }
  else if(gp->getStrnd() == dsu::ePos && gp->getEnd5() == end5 && end5 != -1) {
      const GenePrediction* p = dynamic_cast<const GenePrediction*>(*(dp->begin()));
      assert(p);
      bool isInit = ((p->getScoreIdx() == Initial::getInstance().getScoreIdx()) ||
	(p->getScoreIdx() == Single::getInstance().getScoreIdx())) &&
	p->getEnd5()==end5;
      // make sure not a start
      if( !(isInit && p->getEnd5() == gp->getEnd5()) ) {
	getValRef() = gp->getScore().getVal();
      }
  } else if(gp->getStrnd() == dsu::eNeg && gp->getEnd3() == end3 && end3 != -1) {
      const GenePrediction* p = dynamic_cast<const GenePrediction*>(dp->getLast());
      assert(p);
      bool isInit = ((p->getScoreIdx() == Initial::getInstance().getScoreIdx()) ||
	(p->getScoreIdx() == Single::getInstance().getScoreIdx())) &&
	p->getEnd3()==end3;
      // make sure not start
      if( !(isInit && p->getEnd3() == gp->getEnd3()) ) {
	getValRef() = gp->getScore().getVal();
      }
  }
}

void GenePredEvidenceBndryDon::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  assert(dp);
  const GenePrediction* gp = TypeManage::checkForGenePrediction(pred,end5,end3);
  assert(gp);
  if( gp->getStrnd() != strnd || gp->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
  } else if(gp->getStrnd() == dsu::ePos && gp->getEnd3() == end3 && end3 != -1) {
      const GenePrediction* p = dynamic_cast<const GenePrediction*>(dp->getLast());
      assert(p);
      bool isTerm = ((p->getScoreIdx() == Terminal::getInstance().getScoreIdx()) ||
	(p->getScoreIdx() == Single::getInstance().getScoreIdx())) &&
	p->getEnd3()==end3;
      // make sure not a stop
      if( !(isTerm && p->getEnd3() == gp->getEnd3()) ) {
	getValRef() = gp->getScore().getVal();
      }
    } else if(gp->getStrnd() == dsu::eNeg && gp->getEnd5() == end5 && end5 != -1) {
      const GenePrediction* p = dynamic_cast<const GenePrediction*>(*(dp->begin()));
      assert(p);
      bool isTerm = ((p->getScoreIdx() == Terminal::getInstance().getScoreIdx()) ||
	(p->getScoreIdx() == Single::getInstance().getScoreIdx())) &&
	p->getEnd5()==end5;
      // make sure not a stop
      if( !(isTerm && p->getEnd5() == gp->getEnd5()) ) {
	getValRef() = gp->getScore().getVal();
      }
    }
}

void GenePredEvidenceBndryStart::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd ) {
    return;
  }
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  assert(dp);
  if( pred->getStrnd() == dsu::ePos ) {
    const GenePrediction* p = dynamic_cast<const GenePrediction*>(*(dp->begin()));
    assert(p);
    if( p->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
      return;
    }
    bool isInit = (p->getScoreIdx() == Initial::getInstance().getScoreIdx()) ||
      (p->getScoreIdx() == Single::getInstance().getScoreIdx());
    if( isInit && p->getEnd5() == pred->getEnd5() && p->getStrnd() == dsu::ePos && 
	p->getEnd5() == end5 && end5 != -1)
      getValRef() = p->getScore().getVal();
  } else {
    const GenePrediction* p = dynamic_cast<const GenePrediction*>(dp->getLast());
    assert(p);
    if( p->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
      return;
    }
    bool isInit = (p->getScoreIdx() == Initial::getInstance().getScoreIdx()) ||
      (p->getScoreIdx() == Single::getInstance().getScoreIdx());
    if( isInit && p->getEnd3() == pred->getEnd3() && p->getStrnd() == dsu::eNeg && p->getEnd3() == end3 && end3 != -1) 
      getValRef() = p->getScore().getVal();
  }
}

void GenePredEvidenceBndryStop::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  assert(dp);
  if( pred->getStrnd() == dsu::ePos ) {
    const GenePrediction* p = dynamic_cast<const GenePrediction*>(dp->getLast());
    assert(p);
    if( p->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
      return;
    }
    bool isTerm = (p->getScoreIdx() == Terminal::getInstance().getScoreIdx()) ||
      (p->getScoreIdx() == Single::getInstance().getScoreIdx());
    if( isTerm && p->getEnd3() == pred->getEnd3() && p->getStrnd() == dsu::ePos && p->getEnd3() == end3 && end3 != -1)
      getValRef() = p->getScore().getVal();
  } else {
    const GenePrediction* p = dynamic_cast<const GenePrediction*>(*(dp->begin()));
    assert(p);
    if( p->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
      return;
    }
    bool isTerm = (p->getScoreIdx() == Terminal::getInstance().getScoreIdx()) ||
      (p->getScoreIdx() == Single::getInstance().getScoreIdx());
    if( isTerm && p->getEnd5() == pred->getEnd5() && p->getStrnd() == dsu::eNeg && p->getEnd5() == end5 && end5 != -1) 
      getValRef() = p->getScore().getVal();
  }
}

EvidenceType* GenePredEvidenceIntron::dup() const {
  return new GenePredEvidenceIntron(getStr(),getIdx(),_wght);
}

void GenePredEvidenceIntron::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const GenePrediction* gp = TypeManage::checkForGenePrediction(pred,end5,end3);
  assert(gp);
  if( gp->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
    double val = gp->getScore().getVal();
    getValRef() = val;
  }
}

EvidenceType* GenePredEvidenceOverlapStart::dup() const {
  return new GenePredEvidenceOverlapStart(getStr(),getIdx(),_wght);
}

void GenePredEvidenceOverlapStart::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const GenePrediction* gp = TypeManage::checkForGenePrediction(pred,end5,end3);
  assert(gp);
  // this checking is redundant
  if( (pred->getStrnd()==dsu::ePos && (end5 >= pred->getEnd5() && end5 <= pred->getEnd3())) ||
      (pred->getStrnd()==dsu::eNeg && (end3 >= pred->getEnd5() && end3 <= pred->getEnd3()))) {
  //if( gp->getScoreIdx() == Intron::getInstance().getScoreIdx() ) {
    const double val = gp->getScore().getVal();
    getValRef() = val;
  }
}
