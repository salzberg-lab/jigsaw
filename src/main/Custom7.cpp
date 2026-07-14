//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Custom7.h"
#include "ExonDistr.h"
#include "AnnotationItem.h"
#include "Options.h"
#include "DataDescrip.h"
#include "RegionPrediction.h"
#include <vector>
#include <algorithm>

#include <fstream>

#if __GNUC__ > 3 || ( __GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#include <cmath>
using std::log;
#else
#include <math.h>
#endif

using std::cerr;
using std::endl;
using std::cout;
using std::runtime_error;

#if _RUNTIME_DEBUG >= 1
bool MYDEBUG = false;
#endif
Custom7* Custom7::_instance = NULL;

static AnnotationScore getProb(const RegionPrediction&, EvScores::Bnd_t, dsu::Strand_t, int); 

void Custom7::createTree(const string& rstr) {
  const Options& opt = *Options::getInstance();
  if( opt.isLinComb() ) return;
  string rdir = rstr;
  if( rdir[rdir.size()-1] != '/' ) {
    rdir += '/';
  }
  const int num_vecs = 3;
  if( rstr.size() ) {
    const string astr = rdir + "acc.dt";
    const string astrv = rdir + "acc.con";
    const DataDescrip& descripAcc = *(getFeatMaps()[EvScores::eAcc]);
    _dTree[EvScores::eAcc] = new DTree(astr,astrv,descripAcc,num_vecs);

    const string dstr = rdir + "don.dt";
    const string dstrv = rdir + "don.con";
    const DataDescrip& descripDon = *(getFeatMaps()[EvScores::eDon]); 
    _dTree[EvScores::eDon] = new DTree(dstr,dstrv,descripDon,num_vecs);

    const string cstr = rdir + "code.dt";
    const string cstrv = rdir + "code.con";
    const DataDescrip& descripCod = *(getFeatMaps()[EvScores::eCoding]); 
    _dTree[EvScores::eCoding] = new DTree(cstr,cstrv,descripCod,num_vecs);

    const string sstr = rdir + "start.dt";
    const string sstrv = rdir + "start.con";
    const DataDescrip& descripStart= *(getFeatMaps()[EvScores::eStart]); 
    _dTree[EvScores::eStart] = new DTree(sstr,sstrv,descripStart,num_vecs);
    
    const string pstr = rdir + "stop.dt";
    const string pstrv = rdir + "stop.con";
    const DataDescrip& descripStop= *(getFeatMaps()[EvScores::eStop]); 
    _dTree[EvScores::eStop] = new DTree(pstr,pstrv,descripStop,num_vecs);
    
    if( opt.useIntron() ) {
      const string istr = rdir + "intron.dt";
      const string istrv = rdir + "intron.con";
      const DataDescrip& descripIntron= *(getFeatMaps()[EvScores::eIntron]); 
      _dTree[EvScores::eIntron] = new DTree(istr,istrv,descripIntron,num_vecs);
    }
  }
}

CustomScore::LabelScore Custom7::score(const RegionPrediction& rp, const string& str) const {
  assert(0);
  CustomScore::LabelScore vals(2);
  return vals;
}

#define CBEST

void Custom7::score(vector<RegionPrediction*>& vrp, const string& str) const {
  int middle = (vrp.size()-1)/2;
  // this method considers only physically adjacent vectors
  // when this is on the distance metric will not be used
  if( vrp[middle-1] ) {
    int leftGap = vrp[middle]->getEnd5()-vrp[middle-1]->getEnd3();
    vrp[middle-1] = leftGap > 1 ? NULL : vrp[middle-1]; 
  }
  if( vrp[middle+1] ) {
    int rghtGap = vrp[middle+1]->getEnd5()-vrp[middle]->getEnd3();
    vrp[middle+1] = rghtGap > 1 ? NULL : vrp[middle+1];
  }
  CustomScore::LabelScore vals(2);
  vals[dsu::eNeg] = scoreHelp(vrp,str,dsu::eNeg);
  vals[dsu::ePos] = scoreHelp(vrp,str,dsu::ePos);
  RegionPrediction& rp = *vrp[middle];
  if( !rp.noStrndData(dsu::ePos)) {
    rp.setScores(dsu::ePos,vals[dsu::ePos]);
  }
  if( !rp.noStrndData(dsu::eNeg)) {
    rp.setScores(dsu::eNeg,vals[dsu::eNeg]);
  }
}

EvScores Custom7::scoreHelp(const vector<RegionPrediction*>& rpv, const string& str, dsu::Strand_t strnd) const {
  // we will assume (for now) that the 
  // we are scoring the middle rp of a odd number of vectors
  int middle = (rpv.size()-1)/2;
  if( rpv[middle]->noStrndData(strnd) ) return EvScores();
#if _RUNTIME_DEBUG >= 1
  bool localDebug = false;
  if((rpv[middle]->getEnd5() == 6222)
     && strnd == dsu::ePos ) {
    localDebug = true;
    cout<<"DEBUGON: "<<rpv[middle]->getEnd5()<<" "<<rpv[middle]->getEnd3()<<endl;
    //cout<<*rpv[middle]<<endl;
  } else {
    MYDEBUG = false;
  }
#endif
  double leftGap = 1;
  if( rpv[middle-1] ) {
    leftGap = rpv[middle]->getEnd5()-rpv[middle-1]->getEnd3();
    leftGap = 1 - (1/leftGap);
  } else {
    leftGap = 1;
  }
  double rghtGap = 1;
  if( rpv[middle+1] ) {
    rghtGap = rpv[middle+1]->getEnd5()-rpv[middle]->getEnd3();
    rghtGap = 1 - (1/rghtGap);
  } else {
    rghtGap = 1;
  }
  assert(middle >= 1);
  assert(middle+1 < (signed)rpv.size());
  vector< vector<FeatVec*> > vecofvec(rpv.size());
  for(unsigned i = 0; i < vecofvec.size(); ++i) {
#if _RUNTIME_DEBUG >= 3
    if(localDebug && i == 1)
      MYDEBUG=true;
    else
      MYDEBUG=false;
#endif
    vecofvec[i] = calcVec(rpv[i],strnd);
  }
  EvScores scores;
  const Options& opt = *Options::getInstance();
  const unsigned cNum = opt.useIntron() ? EvScores::eNumBnd : EvScores::eNumBnd-1; 
  bool isLinComb = opt.isLinComb();
  for(unsigned i = 0; i < cNum; ++i) {
    if( !_dTree[i] && !isLinComb ) 
      throw runtime_error("didn't read trees");
    vector<FeatVec*> subvec(vecofvec.size()); 
    for(unsigned k = 0; k < vecofvec.size(); ++k) {
	unsigned idx = k;
	// the vectors are assume to always work 5'->3' (so that +/- vectors are the same)
	if(strnd == dsu::eNeg) {
	  idx = vecofvec.size()-k-1;
	}
	subvec[idx] = vecofvec[k][i];
    }
#if _RUNTIME_DEBUG >= 2
    if(localDebug && (i == EvScores::eCoding) ) {
      MYDEBUG = true;
	cout<<"Classify: ";
	subvec[0]->val_print(cout);
	subvec[1]->val_print(cout);
	subvec[2]->val_print(cout);
	cout<<endl;
    } else 
      MYDEBUG = false;
#endif
    scores[i] = _dTree[i]->eval(subvec);
  }
  if( scores[EvScores::eCoding]==0 ) {
    scores[EvScores::eCoding] = AnnotationScore::badVal();
  }
  for(unsigned i = 0; i < vecofvec.size(); ++i) {
    for(unsigned j = 0; j < cNum; ++j) {
      delete vecofvec[i][j];
    }
  }
  return scores;
}

double Custom7::getProb(const AnnotationItem& item, EvScores::Bnd_t cbnd, dsu::Strand_t strnd) const {
  const RegionPrediction& rp = item.getRegionPrediction();
  double prb = ::getProb(rp,cbnd,strnd,0).getVal();
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

static AnnotationScore 
getRawProb(const RegionPrediction& rp, EvScores::Bnd_t bnd, dsu::Strand_t strnd) {
  //when noncoding start with arbitrary strnd, both will be considered
  strnd = dsu::eEither == strnd && bnd == EvScores::eNonCoding ? dsu::ePos : strnd;
  AnnotationScore score = 1.0;
  for(unsigned i = 0; i < EvScores::eNumBnd; ++i) {
    bvec[i] = false;
  }
  bnd2Pred(bnd,bvec);
  const Options& opt = *Options::getInstance();
  const unsigned cNum = opt.useIntron() ? EvScores::eNumBnd : EvScores::eNumBnd-1; 
  for(unsigned i = 0; i < cNum; ++i) {
    AnnotationScore tmp = rp.getScore(i,strnd);
    if( bvec[i] ) {
      if( tmp == AnnotationScore::badVal() ) return tmp;
      score *= tmp;
    } else {
      score *= (tmp==AnnotationScore::badVal() ? 1 :(1-tmp));
    }
  }
  dsu::Strand_t oppStrnd = (strnd==dsu::eNeg ? dsu::ePos : dsu::eNeg);
  for(unsigned i = 0; i < cNum; ++i) {
    AnnotationScore tmp = rp.getScore(i,oppStrnd);
    score *= (tmp==AnnotationScore::badVal() ? 1 : (1-tmp));
  }
  return score;
}

/*
** The parse is a conjunction of probabilities 
** so to get the most probable sequence we would 
** normally minimized the -log sum , becuase I already
** was doing a maximizing score, I'm taking 1-prob
*/
static AnnotationScore 
getProb(const RegionPrediction& rp, EvScores::Bnd_t bnd, dsu::Strand_t strnd, int exonCnt) {
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
  double logProb = (prob == 0 ? 10 : -log(prob));
  double exonDistr = ExonDistr::getInstance().getProb(exonCnt);
  return logProb+exonDistr;
}
