//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include "Options.h"
#include "DnaStr.h"
#include "RegionPrediction.h"
#include "AnnotationScore.h"
#include "Prediction.h"
#include "genemodels.h"
#include "PlaceHolder.h"
#include "CustomScore.h"
#include "ExonDistr.h"
#include "ParamFile.h"

using std::cerr;
using std::endl;
using std::cout;

#define CUTOFF 4
#if _RUNTIME_DEBUG >= CUTOFF
static bool STOP_FLAG = false;
#endif

static const bool
checkForPrediction5(const Prediction* pred, int end5) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return false;
  list<const Prediction*> lst;
  dp->getRelatedPrediction5(end5,lst);
  return !lst.empty();
}

static const bool
checkForPrediction3(const Prediction* pred, int end3) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return false;
  list<const Prediction*> lst;
  dp->getRelatedPrediction3(end3,lst);
  return !lst.empty();
}

static const string& qlabel = "none";
RegionPrediction::RegionPrediction(const SeqLoc& sl, const string& str) :
  DataSetAbstractProduct(sl,-1,qlabel,AnnotationScore::badVal()), 
  _score(2), _scoreOpp(2),
  _negCnt(0), _posCnt(0), _seq(str) {
  for(int i = 0; i < 2; ++i) {
    _accpBndry[i] = _donrBndry[i] = _hasStart[i] = _hasStop[i] = false;
  }
}

RegionPrediction::~RegionPrediction() {
  _score.clear();
  lnkSeqLoc.clear();
}

extern bool TDEBUG;

void RegionPrediction::storeRelaventData(InputData::const_iterator& search_iter, const InputData::const_iterator& stop) {
  InputData::const_iterator iter = search_iter;
  int farthest_right = -1;
  const bool verbose = Options::getInstance()->isVerbose();
  while( iter != stop ) {
    const PlaceHolder* sl = NULL;
    sl = dynamic_cast<const PlaceHolder*>(*iter);
    if(!sl) {
      cerr<<"problems..."<<endl;
      exit(-1);
    }
    //if( TDEBUG ) {
      //cout<<*sl<<endl;
    //}
    if(sl->getEnd3() > farthest_right)
      farthest_right = sl->getEnd3();

    // completeley to the left of region
    if( sl->getEnd3() < getEnd5() ) {
      ++iter;
      if(sl->getEnd3() >= farthest_right || getEnd5() > farthest_right) {
	search_iter = iter;
        continue;
      } else if( sl->getEnd3() < farthest_right) {
	continue;
      } else {
	break;
      }
    }
    // completeley to the right of region
    if( sl->getEnd5() > getEnd3() )
      break;
    
    if( sl->getEnd5() > getEnd5() ) {
      assert( sl->getEnd5() <= getEnd3() );
      int end5 = sl->getEnd5() - 1;
      setEnd3( end5 );
      doubleCheck();
      break;
    } 
    insert(sl);
    if( verbose && ((length() == 1 && size() > 2000) || size() > 8500 )) {
      //cout<<"Note: ignoring additional overlapping data due to the large size.  "<<size()<<endl;
      break;
    }
    ++iter;
  }
}

bool RegionPrediction::noData() const { 
  return (noStrndData(dsu::eNeg) && noStrndData(dsu::ePos));
} 

static const string strType = "RegionPrediction";
const string& RegionPrediction::getTypeStr() const { 
  return strType;
}

void RegionPrediction::setScores(dsu::Strand_t strnd, const EvScores& es) {
  getScores(strnd) = es;
}

void RegionPrediction::setScoresOpp(dsu::Strand_t strnd, const EvScores& es) {
  _scoreOpp[strnd] = es;
}

/****
 **** BUG HERE, I've got val3 going off as true by coincidence
 ** should get picked up as an intron becuase the frames will miss these 
 ** stop codons should take strand and frame into consideration
 */
bool
RegionPrediction::is_igr(const RegionPrediction& left, const RegionPrediction& rght, dsu::Strand_t strnd,int frame) {
 bool val1 = left.hasStop(dsu::ePos) && rght.hasStart(dsu::ePos) && frame==0 && strnd==dsu::ePos;
 bool val2 = left.hasStart(dsu::eNeg) && rght.hasStart(dsu::ePos) && frame==0 && strnd==dsu::ePos;
 bool val3 = left.hasStop(dsu::ePos) && rght.hasStop(dsu::eNeg) && frame==0 && strnd==dsu::eNeg;
 bool val4 = left.hasStart(dsu::eNeg) && rght.hasStop(dsu::eNeg) && frame==0 && strnd==dsu::eNeg;
 return val1||val2||val3||val4;
}

bool
RegionPrediction::is_intron(const RegionPrediction& left, const RegionPrediction& rght,dsu::Strand_t strnd) {
  bool val1 = left.isDonrBndry(strnd) && rght.isAccpBndry(strnd) && strnd==dsu::ePos;
  bool val2 = left.isAccpBndry(strnd) && rght.isDonrBndry(strnd) && strnd==dsu::eNeg;
  bool val3 = val1 || val2;
  return val3;
}


void RegionPrediction::print(std::ostream& os, const CustomScore& cs) const { 
  os<<getTypeStr()<<" ";
  SeqLoc::print(os);
  os<<" size: "<<size()<<" ";
  //os<<endl;
  os<<":";
  os<<getScores(dsu::ePos);
  os<<getScores(dsu::eNeg);
  os<<getScoresOpp(dsu::ePos);
  os<<getScoresOpp(dsu::eNeg);
  os<<"\nAvailable Evidence Sets: "<<endl;
  vector<FeatVec*> posVec = cs.calcVec(this,dsu::ePos);
  for(unsigned i = 0; i < posVec.size(); ++i) {
    posVec[i]->val_print(os);
    os<<endl;
  }
  vector<FeatVec*> negVec = cs.calcVec(this,dsu::eNeg);
  for(unsigned i = 0; i < negVec.size(); ++i) {
    negVec[i]->val_print(os);
    os<<endl;
  }
  DataStorageType::const_iterator iter = begin();
  while( iter != end() ) {
    const SeqLoc& sl = **iter;
    sl.print(os,"\t");
    os<<endl;
    ++iter;
  }
}

bool RegionPrediction::noStrndData(dsu::Strand_t strnd) const {
#if _RUNTIME_DEBUG >= 1
  assert( strnd != dsu::eEither);
#endif
  bool rVal = false;
  const DataStorageType::const_iterator iter = begin();
  //const DataStorageType::const_iterator stop = begin();
  if( (strnd == dsu::ePos && _posCnt == 1) ||
      (strnd == dsu::eNeg && _negCnt == 1)) {
    rVal = ( (*iter)->length() == 1 && (*iter)->getStrnd() == strnd); 
  } else if((strnd == dsu::ePos && _posCnt == 0) ||
	    (strnd == dsu::eNeg && _negCnt == 0)) {
    rVal = true;
  }
  return rVal;
}


void RegionPrediction::insert(const PlaceHolder*  seqLoc) {
  int end3 = seqLoc->getEnd3(); 
  if( end3 < getEnd3() || !isEnd3Set()) {
    setEnd3(end3);
  }
  if( getEnd5()-getEnd3() == 0 && lnkSeqLoc.empty() ) 
    setEnd3(end3);
  assert(getEnd3() != -1);

  DataStorageType::iterator iter = begin(), save_iter = end(), stop = end();
  while( iter != stop) {
    Prediction& sl = **iter;
    if(seqLoc->getEnd5() < sl.getEnd5() && save_iter == stop ) {
      save_iter = iter;
    }
    DataStorageType::iterator del_iter = iter;
    ++iter;
    if(!sl.isActualOverlap(getEnd5(),getEnd3())) 
      lnkSeqLoc.erase(del_iter);
  }
  dsu::Strand_t strnd = seqLoc->getStrnd();
  assert(strnd != dsu::eEither);
  if(strnd == dsu::ePos)
    ++_posCnt;
  else
    ++_negCnt;
  // annoying const issue!!!
  //#warning "Letting a const issue slip for the moment"
  lnkSeqLoc.insert(save_iter, const_cast<Prediction*>(seqLoc->getWholeShabang()));
}

struct eqstr {
  bool operator()(const char* s1, const char* s2) const { return strcmp(s1, s2) == 0; }
};

struct InitBool {
  InitBool() { 
    for(int i = 0; i < 2; i++) _val[i] = false; 
  }
  bool getVal(dsu::Strand_t strnd) const {return _val[strnd];}
  void setVal(dsu::Strand_t strnd, bool val) {
    _val[strnd] = val;
  }
private:
  bool _val[2];
  dsu::Strand_t _strnd;
};

void RegionPrediction::doubleCheck() {
  DataStorageType::iterator iter = begin(), stop = end();
  while( iter != stop ) {
    Prediction& sl = **iter;
    DataStorageType::iterator del_iter = iter;
    ++iter;
    if(!sl.isActualOverlap(getEnd5(),getEnd3())) 
      lnkSeqLoc.erase(del_iter);
  }
}

void RegionPrediction::checkProteinStops(const string& str, dsu::Strand_t strnd) {
  DnaStr dnastr(str,strnd);
  DataStorageType::const_iterator iter = begin(), stop = end();
  int end5 = getEnd5(), end3 = getEnd3();
  while(iter != stop) {
    const Prediction* pred = *iter;
    if(pred->getStrnd() != strnd) {
      ++iter; continue;
    }
    if(end3+2 < (signed)str.length()) {
      if( dnastr.isStop(end3+1) && strnd == dsu::ePos) {
	if( checkForPrediction3(pred,end3)) {
	  _hasStop[strnd] = true;
	}
      } else if(end5-3 >= 0) {
	if(dnastr.isStop(end5-3) && strnd == dsu::eNeg) {
	  if( checkForPrediction5(pred,end5)) {
	    _hasStop[strnd] = true;
	  }
	}
      }
    }
    ++iter;
  }
}

void RegionPrediction::score(const CustomScore& customScore, const string& str) {
  scoreHelp(str,dsu::ePos);
  scoreHelp(str,dsu::eNeg);
  CustomScore::LabelScore coding_score = customScore.score(*this,str);
  if( customScore.parseOnly() )
    return;
  if( !noStrndData(dsu::ePos)) {
    setScores(dsu::ePos,coding_score[dsu::ePos]);
  }
  if( !noStrndData(dsu::eNeg)) {
    setScores(dsu::eNeg,coding_score[dsu::eNeg]);
  }
}

void RegionPrediction::scoreHelp(const string& str, dsu::Strand_t strnd) {
  if( noStrndData(strnd) ) return;
  int end5 = getEnd5(), end3 = getEnd3();
#if _RUNTIME_DEBUG >= CUTOFF
    if( end5 == 187 ) {
      cout<<*this<<endl;
      STOP_FLAG = true;
    } else {
      STOP_FLAG = false;
    }
#endif
    setAcceptorEvidenceBoundary(str,strnd);
    setDonorEvidenceBoundary(str,strnd);
    checkProteinStops(str,strnd);
    DnaStr dnastr(str,strnd);
    // RIGHT SIDE (+3')
    if(end3 - 2 >= 0) {
      if( dnastr.isStop(end3-2) && strnd == dsu::ePos ) {
	_hasStop[strnd] = true;
      } else if( dnastr.isStart(end3-2) && strnd == dsu::eNeg ) {
	_hasStart[strnd] = true;
      }
    } 
    // LEFT SIDE (+5')
    if( dnastr.isStart(end5) && strnd == dsu::ePos ) {
      _hasStart[strnd] = true;
    } else if(dnastr.isStop(end5) && strnd == dsu::eNeg ) {
      _hasStop[strnd] = true;
    } 
}

static inline unsigned bnd2Idx(EvScores::Bnd_t bnd, dsu::Strand_t strnd) {
  const int dupl = 9;
  int idx;
  switch(strnd) {
  case EvScores::eAcc:
    idx = 0;
    break;
  case EvScores::eDon:
    idx = 1;
    break;
  case EvScores::eCoding:
    idx = 2;
    break;
  case EvScores::eStart:
    idx = 3;
    break;
  case EvScores::eStop:
    idx = 4;
    break;
    //Internal
  case (EvScores::eAccH|EvScores::eDonH|EvScores::eCodingH):
    idx = 5;
    break;
    //Single
  case (EvScores::eStartH|EvScores::eCodingH|EvScores::eStopH):
    idx = 6;
    break;
    //Terminal
  case (EvScores::eStopH | EvScores::eCodingH | EvScores::eAccH):
    idx = 7;
    break;
    //Initial
  case (EvScores::eStartH | EvScores::eCodingH | EvScores::eDonH):
    idx = 8;
    break;
  case EvScores::eNonCoding:
    idx = 18;
    return idx;
    break;
  default:
    assert(0);
  }
  int off;
  if(strnd==dsu::eNeg)
    off=1;
  else
    off=0;
  return off*dupl+idx;
}

bool RegionPrediction::isValid(EvScores::Bnd_t bnd, dsu::Strand_t strnd) const {
  const EvScores& ev = getScores(strnd);
  return ev[bnd].isValid();
}

void RegionPrediction::setAcceptorEvidenceBoundary(const string& dnaStr, dsu::Strand_t strnd) {
  DnaStr str(dnaStr,strnd);
  int end3 = getEnd3(), end5 = getEnd5();
  char a = -1, b = -1;
  if(strnd == dsu::ePos && end5 == 1) {
    _accpBndry[strnd] = false;
  } else if(strnd == dsu::ePos && (end5-3)>=0) {
    b = tolower(str[end5-1]);
    a = tolower(str[end5-2]);
  } else if(strnd==dsu::eNeg && end3 == static_cast<signed>(str.length())) {
    _accpBndry[strnd] = false;
  } else if(strnd==dsu::eNeg && (end3+1) < static_cast<signed>(str.length())) {
    a = tolower(str[end3+2]);
    b = tolower(str[end3+1]);
  }
  const ParamFile& pf = (*ParamFile::getInstance());
  const std::list<string>& dlst = pf.getAcceptorList();
  std::list<string>::const_iterator iter = dlst.begin();
  const std::list<string>::const_iterator stop = dlst.end();
  for(; iter != stop; iter++) {
    const string& seq = *iter;
    assert(seq.length()==2);
    if(a == tolower(seq[0]) && b == tolower(seq[1])) {
      _accpBndry[strnd] = true;
      break;
    }
  }
}

void RegionPrediction::setDonorEvidenceBoundary(const string& dnaStr, dsu::Strand_t strnd) {
  DnaStr str(dnaStr,strnd);
  char a = -1, b = -1;
  int end3 = getEnd3(), end5 = getEnd5();
  if(strnd == dsu::ePos && end3 == static_cast<signed>(str.length())) {
    _donrBndry[strnd] = false;
  } else if(strnd == dsu::ePos && end3+1 < static_cast<signed>(str.length())) {
    a = tolower(str[end3+1]);
    b = tolower(str[end3+2]);
  } else if(strnd==dsu::eNeg && end5 == 1) {
    _donrBndry[strnd] = false;
  } else if(strnd==dsu::eNeg && (end5-3)>=0) {
    a = tolower(str[end5-1]);
    b = tolower(str[end5-2]);
  }
  const ParamFile& pf = (*ParamFile::getInstance());
  const std::list<string>& dlst = pf.getDonorList();
  std::list<string>::const_iterator iter = dlst.begin();
  const std::list<string>::const_iterator stop = dlst.end();
  for(; iter != stop; iter++) {
    const string& seq = *iter;
    assert(seq.length()==2);
    if(a == tolower(seq[0]) && b == tolower(seq[1])) {
      _donrBndry[strnd] = true;
      break;
    }
  }
}
