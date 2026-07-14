//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "HomGeneEvidence.h"
#include "ParamFile.h"
#include "TypeManage.h"
#include "genemodels.h"
#include <stdexcept>

HomGeneEvidence::HomGeneEvidence(const string& str, int idx, double val) 
  : EvidenceType(str,idx,val,0) { }

HomGeneEvidence::~HomGeneEvidence() { }

unsigned HomGeneEvidence::getNum() const { return 1; }
static double getRoundVal(double val) {
  val /= 100.0;
  return val;
}

void HomGeneEvidence::finish() {
}

void HomGeneEvidence::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const MatchPrediction* mp = TypeManage::checkForMatchPrediction(pred,end5,end3);
  assert(mp);
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx()) {
    return;
  }
  double val = getRoundVal(mp->getPcntSim());
  if(val > _val)
    _val = val;
}
void HomGeneEvidence::setVals(const EvidenceType& et) {
  const HomGeneEvidence& het = dynamic_cast<const HomGeneEvidence&>(et);
  _val = het._val;
}

#if 1
double HomGeneEvidence::eval() const {
  assert(0);
  return -1;
}
#endif

EvidenceType* HomGeneEvidence::dup() const {
  return new HomGeneEvidence(getStr(),getIdx(),getVal());
}

#if 1
bool HomGeneEvidence::train_print(std::ostream& os) const {
  os<<std::setprecision(2)<<getVal();
  return true;
}
bool HomGeneEvidence::read(std::istream& is) {
  return is>>_val ? true : false;
}
#endif

HomGeneEvidenceBndry::HomGeneEvidenceBndry(const string& str, int idx, double val, const string& seq) :
  HomGeneEvidence(str,idx,val), _theSeq(seq) { }

HomGeneEvidenceBndryAcc::HomGeneEvidenceBndryAcc(const string& str, int idx, double val, const string& seq) :
  HomGeneEvidenceBndry(str,idx,val,seq) { }

HomGeneEvidenceBndryDon::HomGeneEvidenceBndryDon(const string& str, int idx, 
						 double val, const string& seq) : 
  HomGeneEvidenceBndry(str,idx,val,seq) { }

HomGeneEvidenceBndryStart::HomGeneEvidenceBndryStart(const string& str, int idx, double val,
						     const string& seq) : 
  HomGeneEvidenceBndry(str,idx,val,seq) { }

HomGeneEvidenceBndryStop::HomGeneEvidenceBndryStop(const string& str, int idx, 
						   double val, const string& seq) : 
  HomGeneEvidenceBndry(str,idx,val,seq) { }

void HomGeneEvidenceBndry::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  // this is not used is it?
  throw std::runtime_error("unexpected function call HomGeneEvidenceBndry::update");
}

EvidenceType* HomGeneEvidenceBndry::dup() const {
  return new HomGeneEvidenceBndry(getStr(),getIdx(),getVal(),_theSeq);
}

EvidenceType* HomGeneEvidenceBndryAcc::dup() const {
  return new HomGeneEvidenceBndryAcc(getStr(),getIdx(),getVal(),_theSeq);
}

EvidenceType* HomGeneEvidenceBndryDon::dup() const {
  return new HomGeneEvidenceBndryDon(getStr(),getIdx(),getVal(),_theSeq);
}

EvidenceType* HomGeneEvidenceBndryStart::dup() const {
  return new HomGeneEvidenceBndryStart(getStr(),getIdx(),getVal(),_theSeq);
}

EvidenceType* HomGeneEvidenceBndryStop::dup() const {
  return new HomGeneEvidenceBndryStop(getStr(),getIdx(),getVal(),_theSeq);
}

void HomGeneEvidenceBndryAcc::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  DnaStr str(_theSeq,pred->getStrnd());
  const MatchPrediction* mp = NULL;
  if(pred->getStrnd() == dsu::ePos && end5 != -1 && (end5-3) >= 0 ) {
    mp = TypeManage::checkForMatchPrediction5(pred,end5);
  } else if( pred->getStrnd() == dsu::eNeg && end3 != -1 && (end3+2) < static_cast<signed>(str.length()) ) {
    mp = TypeManage::checkForMatchPrediction3(pred,end3);
  }
  if( !mp ) {
    return;
  }
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx()) {
    return;
  }
  char x1='\0',x2='\0';
  if(pred->getStrnd() == dsu::ePos && mp->getEnd5() == end5 && end5 != -1 && (end5-3) >= 0 ) {
    x1 = tolower(str[end5-2]);
    x2 = tolower(str[end5-1]);
  } else if( pred->getStrnd() == dsu::eNeg && mp->getEnd3() == end3 && end3 != -1 && (end3+1) < static_cast<signed>(str.length()) ) {
    x1 = tolower(str[end3+2]);
    x2 = tolower(str[end3+1]);
  } 
  bool isAcc = false;
  if( x1 != '\0' && x2 != '\0' ) {
    const ParamFile& pf = (*ParamFile::getInstance());
    const std::list<string>& dlst = pf.getAcceptorList();
    std::list<string>::const_iterator iter = dlst.begin();
    const std::list<string>::const_iterator stop = dlst.end();
    for(; iter != stop; iter++) {
      const string& seq = *iter;
      assert(seq.length()==2);
      if(x1 == tolower(seq[0]) && x2 == tolower(seq[1])) {
	isAcc=true;
	break;
      }
    }
  }
  if(isAcc) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
  }
}

void HomGeneEvidenceBndryDon::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const MatchPrediction* mp = NULL;
  DnaStr str(_theSeq,pred->getStrnd());
  if(pred->getStrnd() == dsu::ePos && end3 != -1 && (end3+1) < static_cast<signed>(str.length()) ) {
    mp = TypeManage::checkForMatchPrediction3(pred,end3);
  } else if( pred->getStrnd() == dsu::eNeg && end5 != -1 && (end5-3) >= 0 ) {
    mp = TypeManage::checkForMatchPrediction5(pred,end5);
  }
  if( !mp ) {
    return;
  }
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx()) {
    return;
  }
  char x1 = '\0', x2 = '\0';
  if(pred->getStrnd() == dsu::ePos && mp->getEnd3() == end3 && end3 != -1 && (end3+1) < static_cast<signed>(str.length()) ) {
    x1 = tolower(str[end3+1]);
    x2 = tolower(str[end3+2]);
  } else if( pred->getStrnd() == dsu::eNeg && mp->getEnd5() == end5 && end5 != -1 && (end5-3) >= 0 ) {
    x2 = tolower(str[end5-2]);
    x1 = tolower(str[end5-1]);
  } 
  bool isDon = false;
  if( x1 != '\0' && x2 != '\0' ) {
    const ParamFile& pf = (*ParamFile::getInstance());
    const std::list<string>& dlst = pf.getDonorList();
    std::list<string>::const_iterator iter = dlst.begin();
    const std::list<string>::const_iterator stop = dlst.end();
    for(; iter != stop; iter++) {
      const string& seq = *iter;
      assert(seq.length()==2);
      if(x1 == tolower(seq[0]) && x2 == tolower(seq[1])) {
	isDon = true;
	break;
      }
    }
  }
  if(isDon) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
  }
}

void HomGeneEvidenceBndryStart::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const MatchPrediction* mp = TypeManage::checkForMatchPrediction(pred,end5,end3);
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx()) {
    return;
  }
  DnaStr str(_theSeq,pred->getStrnd());
  bool isStart = false;
  assert(end5 != -1);
  assert(end3 != -1);
  if(pred->getStrnd() == dsu::ePos && mp->getEnd5() == end5) {
    isStart = str.isStart(end5);
  } else if( pred->getStrnd() == dsu::eNeg && mp->getEnd3() == end3 && (end3-3) >= 0 ) {
    isStart = str.isStart(end3-2);
  } 
  if(isStart) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
  }
}

void HomGeneEvidenceBndryStop::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  if( pred->getStrnd() != strnd) {
    return;
  }
  const MatchPrediction* mp = TypeManage::checkForMatchPrediction(pred,end5,end3);
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx()) {
    return;
  }
  DnaStr str(_theSeq,pred->getStrnd());
  bool isStop = false;
  assert(end5 != -1);
  assert(end3 != -1);
  if(pred->getStrnd() == dsu::ePos && mp->getEnd3() == end3 && (end3-3) >= 0) {
    isStop = str.isStop(end3-2);
  } else if( pred->getStrnd() == dsu::eNeg && mp->getEnd5() == end5 ) {
    isStop = str.isStop(end5);
  } 
  if(isStop) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
  }
}

HomGeneEvidenceIntron::HomGeneEvidenceIntron(const string& str, int idx, double val, const string& seq) :
  HomGeneEvidenceBndry(str,idx,val,seq) { }

EvidenceType* HomGeneEvidenceIntron::dup() const {
  HomGeneEvidenceIntron* rVal = new HomGeneEvidenceIntron(getStr(),getIdx(),getVal(),_theSeq);
  return rVal;
}

/***
 *  WARNING!! Does not check to see if the alignments end on consensus splice sites
 *  may or may not be a good thing (probably isn't)
 */

void HomGeneEvidenceIntron::update(const Prediction* pred, int iEnd5, int iEnd3, dsu::Strand_t strnd) {
  const MatchPrediction* mp = TypeManage::checkForMatchPrediction(pred,iEnd5,iEnd3);
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx() && mp->getStrnd()==strnd ) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
  }
}
