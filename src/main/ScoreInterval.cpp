#include "ScoreInterval.h"
#include "ParamFile.h"
#include "AnnotationItem.h"
#include "MatrixIdx.h"
#include "genemodels.h"
#include "CustomScore.h"
#include "TransProb.h"
#include <iostream>
#if __GNUC__ >= 3 && __GNUC_MINOR__ >= 4 || __GNUC__ >= 4
#include <cmath>
using std::log;
#endif

using std::cout;
using std::endl;

using namespace dsu;

#if _RUNTIME_DEBUG >= 2
extern bool MYDEBUG1;
#endif

static AnnotationScore
igrHelpPenalty(int, const AnnotationScore&);
static AnnotationScore
intronHelpPenalty(int, const AnnotationScore&);

static AnnotationScore
internalExonHelpPenalty(int, const AnnotationScore&);

static AnnotationScore
internalExonPenalty(const RegionPrediction&, const RegionPrediction&, const AnnotationScore&);

AnnotationScore 
ScoreInterval::score(const AnnotationItem& lItem, const MatrixIdx& lIdx, 
		     const AnnotationItem& rItem, const MatrixIdx& rIdx, 
		     const int type, const AnnotationScore& inv, const CustomScore& cs) {

  // transition between two states
  const pair<int,dsu::Strand_t> ptype = lItem.getPrevExonType(lIdx);
  double tprob = 0;
  const TransProb& tp = TransProb::getInstance();
  if( &lItem != &rItem ) {
    dsu::Strand_t strnd = rIdx.getStrnd();
    tprob = tp.getProb(ptype.first,ptype.second,type,strnd);
  } 
  // initial state
  else {
    // not implemented yet
    //tprob = tp.getProb(type,strnd);
    tprob = 1;
  }
  const int len = (signed)lItem.getSeq(dsu::ePos).length();
  const AnnotationScore storedScore = lItem.getBest(lIdx).getScore();
  AnnotationScore invScore = inv;
  EvScores::Bnd_t bnd = EvScores::eErr, lbnd = EvScores::eErr, rbnd = EvScores::eErr;
  if( type == Single::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eStartH|EvScores::eStopH);
    } else if( lIdx.getStrnd()==dsu::ePos ) {
      assert(rIdx.getStrnd()==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStartH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStopH : EvScores::eCodingH;
    } else {
      assert(lIdx.getStrnd()==dsu::eNeg);
      assert(rIdx.getStrnd()==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStopH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStartH : EvScores::eCodingH;
    }
  } else if( type == Terminal::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eAccH|EvScores::eStopH);
    } else if( lIdx.getStrnd()==dsu::ePos ) {
      assert(rIdx.getStrnd()==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eAccH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStopH : EvScores::eCodingH;
    } else {
      assert(lIdx.getStrnd()==dsu::eNeg);
      assert(rIdx.getStrnd()==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStopH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eAccH : EvScores::eCodingH;
    }
  } else if( type == Initial::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eStartH|EvScores::eDonH);
    } else if( lIdx.getStrnd()==dsu::ePos ) {
      assert(rIdx.getStrnd()==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStartH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eDonH : EvScores::eCodingH;
    } else {
      assert(lIdx.getStrnd()==dsu::eNeg);
      assert(rIdx.getStrnd()==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eDonH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStartH : EvScores::eCodingH;
    }
  } else if( type == Internal::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eAccH|EvScores::eDonH);
    } else if( lIdx.getStrnd()==dsu::ePos ) {
      assert(rIdx.getStrnd()==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eAccH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eDonH : EvScores::eCodingH;
    } else {
      assert(lIdx.getStrnd()==dsu::eNeg);
      assert(rIdx.getStrnd()==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eDonH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eAccH : EvScores::eCodingH;
    }
  }
  AnnotationScore currs = 0;
  if( bnd != EvScores::eErr ) {
    currs = cs.getProb(lItem,bnd,lIdx.getStrnd()); 
    // compute prob that intron doesn't occur
  } else if( lbnd != EvScores::eErr && lbnd != EvScores::eErr ) {
    const AnnotationScore currs1 = cs.getProb(lItem,lbnd,lIdx.getStrnd()); 
    const AnnotationScore currs2 = cs.getProb(rItem,rbnd,rIdx.getStrnd()); 
    if(currs1 != AnnotationScore::badVal() && currs2 != AnnotationScore::badVal()) {
      currs = currs1 + currs2;
    }
  } else if( type == Intron::getInstance().getScoreIdx() ) {
    assert(lIdx.getSide()==eRght && rIdx.getSide()==eLeft);
    const AnnotationScore iscore = inv;
    invScore = intronPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
  } else if( type == Igr::getInstance().getScoreIdx() && lIdx.getStrnd() == rIdx.getStrnd()) {
    const AnnotationScore iscore = inv;
    invScore = igrPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
  }
  if( type == Internal::getInstance().getScoreIdx() ) {
    //const AnnotationScore iscore = inv+currs;
    //currs = internalExonPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
    const AnnotationScore iscore = inv;
    invScore = internalExonPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
  }
  assert(storedScore != AnnotationScore::badVal());
  assert(tprob > 0);
  const double transScore = -log(tprob);
  return currs+invScore+storedScore+transScore;
}

AnnotationScore 
ScoreInterval::scorePartial1ToRight(const AnnotationItem& rItem, const MatrixIdx& rIdx,
				    const int type, const AnnotationScore& inv, const CustomScore& cs) {
  EvScores::Bnd_t rbnd = EvScores::eErr;
  assert(type != Single::getInstance().getScoreIdx());
  if( type == Terminal::getInstance().getScoreIdx() && rIdx.getStrnd()==dsu::ePos ) {
    rbnd = EvScores::eStopH;
  } else if( type == Initial::getInstance().getScoreIdx() && rIdx.getStrnd() == dsu::eNeg ) {
    rbnd = EvScores::eStartH;
  } else if( type == Internal::getInstance().getScoreIdx() && rIdx.getStrnd() == dsu::ePos ) {
    rbnd = EvScores::eDonH;
  } else if( type == Internal::getInstance().getScoreIdx() && rIdx.getStrnd() == dsu::eNeg ) {
    rbnd = EvScores::eAccH;
  }
  assert(rbnd != EvScores::eErr);
  AnnotationScore currs = cs.getProb(rItem,rbnd,rIdx.getStrnd()); 
  currs += inv;
  return currs;
}

AnnotationScore 
ScoreInterval::scorePartialToEnd(const AnnotationItem& rItem, const MatrixIdx& rIdx,
				 const int type, const AnnotationScore& inv) {
  bool intronMod = false;
  if( type == Terminal::getInstance().getScoreIdx() && rIdx.getStrnd()==dsu::eNeg ) {
    intronMod=true;
  } else if( type == Initial::getInstance().getScoreIdx() && rIdx.getStrnd() == dsu::ePos ) {
    intronMod=true;
  } else if( type == Internal::getInstance().getScoreIdx() && (rIdx.getStrnd() == dsu::ePos || rIdx.getStrnd() == dsu::eNeg) ) {
    intronMod=true;
  }
  AnnotationScore nonCode = inv;
  if(intronMod ) {
    const int rghtCoord = rItem.getSeq(dsu::ePos).length();
    const int leftCoord = rItem.getEnd3();
    nonCode = intronPenalty(leftCoord,rghtCoord,inv);
  }
  const AnnotationScore val = rItem.getBest(rIdx).getScore();
  const AnnotationScore currs = nonCode+val;
  return currs;
}


AnnotationScore 
ScoreInterval::igrPenalty(const RegionPrediction& left, const RegionPrediction& rght, const AnnotationScore& score) {
  int dist = (rght.getEnd5()-1)-(left.getEnd3()+1);
  dist += 1;
  return igrHelpPenalty(dist,score);
}


AnnotationScore 
internalExonPenalty(const RegionPrediction& left, const RegionPrediction& rght, const AnnotationScore& score) {
  int dist = rght.getEnd3()-left.getEnd5();
  dist += 1;
  assert( dist >= 1);
  return internalExonHelpPenalty(dist,score);
}

AnnotationScore 
ScoreInterval::intronPenalty(const RegionPrediction& left, const RegionPrediction& rght, const AnnotationScore& score) {
  int dist = (rght.getEnd5()-1)-(left.getEnd3()+1);
  dist += 1;
  return intronHelpPenalty(dist,score);
}

AnnotationScore 
ScoreInterval::intronPenalty(const int left,  const int rght, const AnnotationScore& score) {
  int dist = (rght-1)-(left+1);
  dist += 1;
  return intronHelpPenalty(dist,score);
}

AnnotationScore
intronHelpPenalty(int dist, const AnnotationScore& score) {
  const ParamFile& pf = *(ParamFile::getInstance());
  AnnotationScore val = score;
  // default threshold = 0
  int threshold = pf.getIntronLengthPenalty();
  // default factor = 0.0
  double factor = pf.getIntronPenalty();
  if( dist > threshold && threshold > 0 && factor > 0) {
    val = (factor*(score+1))+score;
  }
  return val;
}


AnnotationScore
igrHelpPenalty(int dist, const AnnotationScore& score) {
  const ParamFile& pf = *(ParamFile::getInstance());
  AnnotationScore val = score;
  // default threshold = 0
  int threshold = pf.getIgrLengthPenalty();
  // default factor = 0.0
  double factor = pf.getIgrPenalty();
  if( dist < threshold && threshold > 0 && factor > 0) {
    val = (factor*(score+1))+score;
  }
  return val;
}

AnnotationScore
internalExonHelpPenalty(int dist, const AnnotationScore& score) {
  const ParamFile& pf = *(ParamFile::getInstance());
  AnnotationScore val = score;
  // default threshold = 0
  int threshold = pf.getInternalExonLengthPenalty();
  // default factor = 0.0
  double factor = pf.getInternalExonPenalty();
  if( dist < threshold && threshold > 0 && factor > 0) {
    val = (factor*(score+1))+score;
  }
  return val;
}


AnnotationScore 
ScoreInterval::scoreOne(const AnnotationItem& lItem, const AnnotationItem& rItem, const int type, dsu::Strand_t strnd, const AnnotationScore& inv, const CustomScore& cs) {
  const int len = (signed)lItem.getSeq(dsu::ePos).length();
  AnnotationScore invScore = inv;
  EvScores::Bnd_t bnd = EvScores::eErr, lbnd = EvScores::eErr, rbnd = EvScores::eErr;
  if( type == Single::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eStartH|EvScores::eStopH);
    } else if( strnd==dsu::ePos ) {
      assert(strnd==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStartH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStopH : EvScores::eCodingH;
    } else {
      assert(strnd==dsu::eNeg);
      assert(strnd==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStopH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStartH : EvScores::eCodingH;
    }
  } else if( type == Terminal::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eAccH|EvScores::eStopH);
    } else if( strnd==dsu::ePos ) {
      assert(strnd==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eAccH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStopH : EvScores::eCodingH;
    } else {
      assert(strnd==dsu::eNeg);
      assert(strnd==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStopH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eAccH : EvScores::eCodingH;
    }
  } else if( type == Initial::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eStartH|EvScores::eDonH);
    } else if( strnd==dsu::ePos ) {
      assert(strnd==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eStartH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eDonH : EvScores::eCodingH;
    } else {
      assert(strnd==dsu::eNeg);
      assert(strnd==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eDonH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eStartH : EvScores::eCodingH;
    }
  } else if( type == Internal::getInstance().getScoreIdx() ) {
    if( &lItem == &rItem ) {
      bnd = (EvScores::Bnd_t)(EvScores::eAccH|EvScores::eDonH);
    } else if( strnd==dsu::ePos ) {
      assert(strnd==dsu::ePos);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eAccH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eDonH : EvScores::eCodingH;
    } else {
      assert(strnd==dsu::eNeg);
      assert(strnd==dsu::eNeg);
      lbnd = lItem.getEnd5() > 1 ? EvScores::eDonH : EvScores::eCodingH;
      rbnd = rItem.getEnd3() < len ? EvScores::eAccH : EvScores::eCodingH;
    }
  }
  AnnotationScore currs = 0;
  if( bnd != EvScores::eErr ) {
    currs = cs.getProb(lItem,bnd,strnd); 
    // compute prob that intron doesn't occur
  } else if( lbnd != EvScores::eErr && lbnd != EvScores::eErr ) {
    const AnnotationScore currs1 = cs.getProb(lItem,lbnd,strnd); 
    const AnnotationScore currs2 = cs.getProb(rItem,rbnd,strnd); 
    if(currs1 != AnnotationScore::badVal() && currs2 != AnnotationScore::badVal()) {
      currs = currs1 + currs2;
    }
  } else if( type == Intron::getInstance().getScoreIdx() ) {
    //assert(lIdx.getSide()==eRght && rIdx.getSide()==eLeft);
    const AnnotationScore iscore = inv; // compute
    invScore = intronPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
  } else if( type == Igr::getInstance().getScoreIdx() ) {
    const AnnotationScore iscore = inv; // compute
    invScore = igrPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
  }
  if( type == Internal::getInstance().getScoreIdx() ) {
    const AnnotationScore iscore = inv+currs;
    currs = internalExonPenalty(lItem.getRegionPrediction(),rItem.getRegionPrediction(),iscore);
    //cout<<"humm: "<<iscore<<" "<<currs<<endl;
  }
  return currs+invScore;
}
