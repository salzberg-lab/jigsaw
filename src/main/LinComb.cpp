//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "LinComb.h"
#include "ExonDistr.h"
#include "AnnotationItem.h"
#include "Options.h"
#include "DataDescrip.h"
#include "RegionPrediction.h"
#include <vector>
#include <algorithm>
#include <fstream>
#if ___GNUC__ > 3 || ( __GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#include <cmath>
using std::log;
#else
#include <math.h>
#endif

using std::cerr;
using std::endl;
using std::cout;
using std::runtime_error;

LinComb* LinComb::_instance = NULL;

//static AnnotationScore getProb(const RegionPrediction&, EvScores::Bnd_t, dsu::Strand_t, int); 

CustomScore::LabelScore LinComb::score(const RegionPrediction& rp, const string& str) const {
  assert(0);
  CustomScore::LabelScore vals(2);
  return vals;
}

#define CBEST

void LinComb::score(vector<RegionPrediction*>& vrp, const string& str) const {
  // this method considers only physically adjacent vectors
  // when this is on the distance metric will not be used
  CustomScore::LabelScore vals(2);
  vals[dsu::eNeg] = scoreHelp(vrp,str,dsu::eNeg);
  vals[dsu::ePos] = scoreHelp(vrp,str,dsu::ePos);
  CustomScore::LabelScore oppvals(2);
  oppvals[dsu::eNeg] = scoreHelpOpp(vrp,str,dsu::eNeg);
  oppvals[dsu::ePos] = scoreHelpOpp(vrp,str,dsu::ePos);
  assert(vrp[1]);
  RegionPrediction& rp = *vrp[1];
  if( !rp.noStrndData(dsu::ePos)) {
    rp.setScores(dsu::ePos,vals[dsu::ePos]);
    rp.setScoresOpp(dsu::ePos,oppvals[dsu::ePos]);
  }
  if( !rp.noStrndData(dsu::eNeg)) {
    rp.setScores(dsu::eNeg,vals[dsu::eNeg]);
    rp.setScoresOpp(dsu::eNeg,oppvals[dsu::eNeg]);
  }
}

EvScores LinComb::scoreHelp(const vector<RegionPrediction*>& rpv, const string& str, dsu::Strand_t strnd) const {
  // we will assume (for now) that the 
  // we are scoring the middle rp of a odd number of vectors
  int middle = (rpv.size()-1)/2;
  if( rpv[middle]->noStrndData(strnd) ) return EvScores();
  assert(middle >= 1);
  assert(middle+1 < (signed)rpv.size());
  vector<FeatVec*> fv = calcVec(rpv[middle],strnd);
  EvScores scores;
  const Options& opt = *Options::getInstance();
  const unsigned cNum = opt.useIntron() ? EvScores::eNumBnd : EvScores::eNumBnd-1; 
  for(unsigned i = 0; i < cNum; ++i) {
    scores[i] = fv[i]->voteEval();
  }
  if( scores[EvScores::eCoding]==0 ) {
    scores[EvScores::eCoding] = AnnotationScore::badVal();
  }
  return scores;
}

EvScores LinComb::scoreHelpOpp(const vector<RegionPrediction*>& rpv, const string& str, dsu::Strand_t strnd) const {
  // we will assume (for now) that the 
  // we are scoring the middle rp of a odd number of vectors
  int middle = (rpv.size()-1)/2;
  if( rpv[middle]->noStrndData(strnd) ) return EvScores();
  assert(middle >= 1);
  assert(middle+1 < (signed)rpv.size());
  const RegionPrediction* currRp = rpv[middle];
  vector<FeatVec*> fv = calcVec(currRp,strnd);
  EvScores scores;
  const Options& opt = *Options::getInstance();
  const unsigned cNum = opt.useIntron() ? EvScores::eNumBnd : EvScores::eNumBnd-1; 
  for(unsigned i = 0; i < cNum; ++i) {
    FeatVec* fl = fv[i];
    scores[i] = fl->voteEvalOpp();
  }
  if( scores[EvScores::eCoding]==0 ) {
    scores[EvScores::eCoding] = AnnotationScore::badVal();
  }
  return scores;
}

double LinComb::getProb(const AnnotationItem& item, EvScores::Bnd_t cbnd, dsu::Strand_t strnd) const {
  const RegionPrediction& rp = item.getRegionPrediction();
  double prb = getProbHelp(rp,cbnd,strnd,0).getVal();
  return prb;
}

static inline void 
bnd2Pred(EvScores::Bnd_t inBnd, vector<bool>& bvec) {
  assert(inBnd != EvScores::eNumBnd);
  assert(bvec.size()==EvScores::eNumBnd);
  //initialize to false by caller
  if( inBnd < EvScores::eNumBnd ) {
    // if only one boundary conditoin
    // was reported, than neccessarily
    // the remaing scores would be coding
    bvec[inBnd] = true;
    bvec[EvScores::eCoding] = (inBnd != EvScores::eIntron);
  } else if(inBnd != EvScores::eNonCoding) {
    assert(inBnd != EvScores::eIntron);
    bvec[EvScores::eStart] = (inBnd & EvScores::eStartH) ? true : false;
    bvec[EvScores::eStop] = (inBnd & EvScores::eStopH) ? true : false;
    bvec[EvScores::eAcc] = (inBnd & EvScores::eAccH) ? true : false;
    bvec[EvScores::eDon] = (inBnd & EvScores::eDonH) ? true : false;
    bvec[EvScores::eCoding] = true;
    bvec[EvScores::eIntron] = false;
  }
}
  /* 
  ** P(l|Vi-1,Vi,Vi+1) 
  **
  ** Let l be a subset of L = {eAcc(a),eDon(d),eCoding(c),eStart(r),eStop(t),eNonCoding(n)}
  ** The labeling l for this sequence interval can be 
  ** any legal subset of L.  Where the rules are defined by the
  ** standard definitions of gene models.  This results in 11 
  ** possible combinations per strand. Each interval always predicts
  ** coding. Therefore, boundary conditions are always
  ** considered in conjunction with an associated coding interval, this
  ** helps compress the number intervals the combiner has to look at.
  ** So if there is no evidence of coding, boundary conditions are currently
  ** ignored.
  **
  ** V is a set of vectors numerically encoding the evidence for the sequence interval
  ** Each element of L is considered to occur independantly of the other elements
  ** of L, a simplifying assumption used to calculated generalized probabilities.
  ** Becuase each element of L is indpendent, it is represented with it's own vector encoding.
  ** A vector for a given element of L is denoted below as V(type) e.g. V(eAcc).
  **
  ** The probability of any given Labeling for interval i, given Vi-1,Vi,Vi+1 is 
  ** the product of the probability of each component of labeling l occuring * the product of the probability
  ** of the remaining lables not occuring.  This includes the evidence for both strands.
  **
  ** Example: P( a ^ c | Vi-1,Vi, Vi+1 ) = P(eAcc|V(a)i-1,V(a)i,V(a)i+1)*
  **                                       P(c|V(c)i-1,V(c)i,V(c)i+1)*
  **                                       1-P(r|V(r)i-1,V(r)i,V(r)i+1)*
  **                                       1-P(t|V(t)i-1,V(t)i,V(t)i+1)*
  **                                       1-P(d|V(d)i-1,V(d)i,V(d)i+1)*
  **                                       1-P(n|V(n)i-1,V(n)i,V(n)i+1)
  **
  ** The assumption to this point is that the position dependant evidence is not a predictor
  ** of non-coding sequence.  Therefore there is no adequet numerical vector encoding representing
  ** non coding, although, it's quite possible that the same representation for coding intervals
  ** is a reasonable representation of non-coding intervals.  If this were the case then it's values
  ** would still just be 1-the coding probilities. Becuase there is no obvious, unique representation
  ** for non coding sequence based on the input evidence alone, the
  ** Non coding intervals must be defined in terms of the other types of evidence. Therefore V does
  ** not include a non coding vector.
  ** P(N|Vi-1,Vi,Vi+1) = 1-P(a)*1-P(d)*1-P(c)*1-P(r)*1-P(t) (vectors not included for shorthand)
  ** Becuase P(N) is dependant on the other 5 labels, it need not be included in the calculation probability
  ** calculatoin due to it's redundancy.
  */


// profiler is telling me to optimize this function
static vector<bool>   bvec(EvScores::eNumBnd,false);

AnnotationScore 
LinComb::getRawProb(const RegionPrediction& rp, EvScores::Bnd_t bnd, dsu::Strand_t strnd) const {
  //when noncoding start with arbitrary strnd, both will be considered
  strnd = dsu::eEither == strnd && bnd == EvScores::eNonCoding ? dsu::ePos : strnd;
  AnnotationScore score = 0.0;
  for(unsigned i = 0; i < EvScores::eNumBnd; ++i) {
    bvec[i] = false;
  }
  bnd2Pred(bnd,bvec);
  const Options& opt = *Options::getInstance();
  const unsigned cNum = opt.useIntron() ? EvScores::eNumBnd : EvScores::eNumBnd-1; 
  const vector<DataDescrip*>& ddVec = getFeatMaps();
  for(unsigned i = 0; i < cNum; ++i) {
    if( bvec[i] ) {
      const AnnotationScore tmp = rp.getScore(i,strnd);
      if( tmp == AnnotationScore::badVal() ) return tmp;
      score += tmp;
    } else {
      const double maxVal = ddVec[i]->getTotWght();
      const AnnotationScore tmp = rp.getScoreOpp(i,strnd);
      score += (tmp==AnnotationScore::badVal() ? (AnnotationScore)maxVal :tmp);
    }
  }
  dsu::Strand_t oppStrnd = (strnd==dsu::eNeg ? dsu::ePos : dsu::eNeg);
  for(unsigned i = 0; i < cNum; ++i) {
    const AnnotationScore tmp = rp.getScoreOpp(i,oppStrnd);
    const double maxVal = ddVec[i]->getTotWght();
    score += (tmp==AnnotationScore::badVal() ? maxVal :tmp);
  }
  return score;
}

/*
** The parse is a conjunction of probabilities 
** so to get the most probable sequence we would 
** normally minimized the -log sum , becuase I already
** was doing a maximizing score, I'm taking 1-prob
*/
AnnotationScore 
LinComb::getProbHelp(const RegionPrediction& rp, EvScores::Bnd_t bnd, dsu::Strand_t strnd, int exonCnt) const {
  assert(exonCnt >= 0);
  AnnotationScore val = getRawProb(rp,bnd,strnd);
  double prob = -1;
  // this occurs when we're looking at the opposite
  // strand, which has 0 evidence
  // in these cases assume 0 propbability
  if( val == AnnotationScore::badVal() ) {
    prob = 0;
  } else {
    prob = val.getVal();
  }
  //double logProb = (prob == 0 ? 10 : -log(prob));
  //double exonDistr = ExonDistr::getInstance().getProb(exonCnt);
  //const int len = rp.getEnd5()-rp.getEnd3()+1;
  return prob;
}
