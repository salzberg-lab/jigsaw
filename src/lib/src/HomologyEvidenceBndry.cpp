//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "HomologyEvidence.h"
#include "TypeManage.h"
#include "ParamFile.h"
#include "genemodels.h"
#include <stdexcept>

  HomologyEvidence::HomologyEvidence(const string& str, int idx, double maxWght, double numWght, double norm) 
    : EvidenceType(str,idx,0,maxWght),
    _sumVal(0), _numWght(numWght), _norm(norm) { }

  HomologyEvidence::~HomologyEvidence() { }

  unsigned HomologyEvidence::getNum() const { return 2; }
  static double getRoundVal(double val) {
    val /= 100.0;
    return val;
  }
  inline double logN(double base, double val) {
    return log10(val)/log10(base);
  }

  void HomologyEvidence::updateSum(double val) {
    _sumVal += val;
  }
void HomologyEvidence::finish() {
  _sumVal /= 100;
  if( _sumVal > 1 )
    _sumVal=1;
  assert(_sumVal >= 0);
  
}

  void HomologyEvidence::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
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
    updateSum(val);
  }
  void HomologyEvidence::setVals(const EvidenceType& et) {
    const HomologyEvidence& het = dynamic_cast<const HomologyEvidence&>(et);
    _val = het._val;
    _sumVal = het._sumVal;
  }

//double HomologyEvidence::eval() const {
    ////return (_maxVal*_maxWght)+(_sumVal*_numWght);
    //return (_maxVal*_wght);
  //}

  EvidenceType* HomologyEvidence::dup() const {
    return new HomologyEvidence(getStr(),getIdx(),_wght,_numWght,_norm);
  }

  bool HomologyEvidence::train_print(std::ostream& os) const {
    os<<std::setprecision(2)<<getVal(0)<<" "<<std::setprecision(2)<<getVal(1);
    return true;
  }
  bool HomologyEvidence::read(std::istream& is) {
    return is>>_val>>_sumVal ? true : false;
  }
  double HomologyEvidence::getVal(int idx) const {
    switch(idx) {
    case 0:
      return _val;
    case 1:
      return _sumVal;
    default:
      std::cerr<<"huh? getval homology.."<<std::endl;
      return _val;
    }
  }

HomologyEvidenceBndry::HomologyEvidenceBndry(const string& str, int idx, 
					     double maxWght, double numWght, double norm, const string& seq) : 
	    HomologyEvidence(str,idx,maxWght,numWght,norm), _theSeq(seq) { }

HomologyEvidenceBndryAcc::HomologyEvidenceBndryAcc(const string& str, int idx, 
					     double maxWght, double numWght, double norm, const string& seq) : 
	    HomologyEvidenceBndry(str,idx,maxWght,numWght,norm,seq) { }

HomologyEvidenceBndryDon::HomologyEvidenceBndryDon(const string& str, int idx, 
					     double maxWght, double numWght, double norm, const string& seq) : 
	    HomologyEvidenceBndry(str,idx,maxWght,numWght,norm,seq) { }

HomologyEvidenceBndryStart::HomologyEvidenceBndryStart(const string& str, int idx, 
					     double maxWght, double numWght, double norm, const string& seq) : 
	    HomologyEvidenceBndry(str,idx,maxWght,numWght,norm,seq) { }

HomologyEvidenceBndryStop::HomologyEvidenceBndryStop(const string& str, int idx, 
					     double maxWght, double numWght, double norm, const string& seq) : 
	    HomologyEvidenceBndry(str,idx,maxWght,numWght,norm,seq) { }

void HomologyEvidenceBndry::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  // this is not used is it?
  throw std::runtime_error("unexpected function call HomologyEvidenceBndry::update");
}

EvidenceType* HomologyEvidenceBndry::dup() const {
  return new HomologyEvidenceBndry(getStr(),getIdx(),_wght,_numWght,_norm,_theSeq);
}

EvidenceType* HomologyEvidenceBndryAcc::dup() const {
  return new HomologyEvidenceBndryAcc(getStr(),getIdx(),_wght,_numWght,_norm,_theSeq);
}

EvidenceType* HomologyEvidenceBndryDon::dup() const {
  return new HomologyEvidenceBndryDon(getStr(),getIdx(),_wght,_numWght,_norm,_theSeq);
}

EvidenceType* HomologyEvidenceBndryStart::dup() const {
  return new HomologyEvidenceBndryStart(getStr(),getIdx(),_wght,_numWght,_norm,_theSeq);
}

EvidenceType* HomologyEvidenceBndryStop::dup() const {
  return new HomologyEvidenceBndryStop(getStr(),getIdx(),_wght,_numWght,_norm,_theSeq);
}

void HomologyEvidenceBndryAcc::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
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
  //const MatchPrediction* mp = TypeManage::checkForMatchPrediction(pred,end5,end3);
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
    updateSum(val);
  }
}

void HomologyEvidenceBndryDon::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
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
  char x1='\0',x2='\0';
  if(pred->getStrnd() == dsu::ePos && mp->getEnd3() == end3 && end3 != -1 && (end3+1) < static_cast<signed>(str.length()) ) {
    x1 = tolower(str[end3+1]);
    x2 = tolower(str[end3+2]);
  } else if( pred->getStrnd() == dsu::eNeg && mp->getEnd5() == end5 && end5 != -1 && (end5-3) >= 0 ) {
    x1 = tolower(str[end5-1]);
    x2 = tolower(str[end5-2]);
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
	isDon=true;
	break;
      }
    }
  }
  if(isDon) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
    updateSum(val);
  }
}

void HomologyEvidenceBndryStart::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
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
    updateSum(val);
  }
}

void HomologyEvidenceBndryStop::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
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
    updateSum(val);
  }
}

HomologyEvidenceIntron::HomologyEvidenceIntron(const string& str, int idx, 
					     double maxWght, double numWght, double norm, const string& seq) : 
	    HomologyEvidenceBndry(str,idx,maxWght,numWght,norm,seq), _intronPred(0), _nonIntronPred(0) { }

EvidenceType* HomologyEvidenceIntron::dup() const {
  HomologyEvidenceIntron* rVal = new HomologyEvidenceIntron(getStr(),getIdx(),_wght,_numWght,_norm,_theSeq);
  rVal->_intronPred = _intronPred;
  rVal->_nonIntronPred = _nonIntronPred;
  return rVal;
}

/***
 *  WARNING!! Does not check to see if the alignments end on consensus splice sites
 *  may or may not be a good thing (probably isn't)
 */

void HomologyEvidenceIntron::update(const Prediction* pred, int end5, int end3, dsu::Strand_t strnd) {
  const MatchPrediction* mp = TypeManage::checkForMatchPrediction(pred,end5,end3);
  assert(mp);
  if( mp->getScoreIdx() == Intron::getInstance().getScoreIdx() && mp->getStrnd()==strnd ) {
    double val = getRoundVal(mp->getPcntSim());
    if(val > _val)
      _val = val;
    ++_intronPred;
  } else {
    ++_nonIntronPred;
  }
}

void HomologyEvidenceIntron::finish() {
  const double sum = (double)(_intronPred+_nonIntronPred);
  _sumVal = sum > 0 ? (double)_intronPred/(double)(_intronPred+_nonIntronPred) : 0;
}
