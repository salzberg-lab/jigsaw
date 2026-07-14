//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "AnnotationItem.h"
#include "ScoreInterval.h"
#include "SingleExon.h"
#include "InitialExon.h"
#include "IntronInterval.h"
#include "IgrInterval.h"
#include "InternalExon.h"
#include "TerminalExon.h"
#include "Compare.h"
#include "Options.h"
#include "RegionPrediction.h"
#include "genemodels.h"
#include "AnnotationScore.h" 
#include "DnaStr.h"
#include "CustomScore.h"
#include "ExonDistr.h"
#include <assert.h>

using std::cout;
using std::endl;
using std::cerr;
using namespace dsu;

#if _RUNTIME_DEBUG >= 2
extern bool DEBUGON;
#endif

/** speed things up a bit... **/
static char justDontAskAboutThisOptimizationHackShit[3];
#if _RUNTIME_DEBUG >= 2
static bool noStops(const DnaStr& seq, int start, int stop, bool verbose) {
#else
static bool noStops(const DnaStr& seq, int start, int stop) {
#endif 
  for(int i = start; i <= stop-2; i+=3) { 
#if _RUNTIME_DEBUG >= 2
    if(verbose)
      seq.printCodon(i);
#endif
    if(seq.isStop(i)) {
      return false;
    }
  } 
  return true;
} 

AnnotationItem::AnnotationItem(const vector<RegionPrediction*>& rpVec, const DnaStr& posSeq, const DnaStr& negSeq) : 
  _seq(*rpVec[1]), _posSeq(posSeq), _negSeq(negSeq), _strnd(dsu::eEither), _best(2*2*6), _left(rpVec[0]), _rght(rpVec[2]) {
}

MatrixItemContents&
AnnotationItem::getBest(const MatrixIdx& mi) {
  return _best[mi.getIntIdx()];
}

const MatrixItemContents&
  AnnotationItem::getBest(const MatrixIdx& mi) const {
  return _best[mi.getIntIdx()];
}

pair<AnnotationScore,MatrixItemContents>
AnnotationItem::endScore(const CustomScore& cs, const AnnotationScore& inScore, bool allowPartialGene) const {
  const Options& opt = *Options::getInstance();
  const Compare* cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    cmp = GtCompare::getInstance();
  AnnotationScore best_score = AnnotationScore::badVal();
  MatrixItemContents tbk;
  dsu::side_t side = eRght;
  for(int bnd= 0; bnd < 2; ++bnd) {
    dsu::gbnd_t gb = (dsu::gbnd_t)bnd;
    if( (allowPartialGene && gb == dsu::eMdlOfGene) || gb == dsu::eEndOfGene ) {
      for(int idx = 0; idx < 2; ++idx) {
	dsu::Strand_t strnd = (dsu::Strand_t)idx;
	for(unsigned frame = 0; frame < 3; ++frame) {
	  MatrixIdx mIdx(frame,strnd,side,gb);
	  if(getBest(mIdx).getScore().isValid()) {
	    const int type = getBest(mIdx).getType();
	    const AnnotationScore val = ScoreInterval::scorePartialToEnd(*this,mIdx,type,inScore);
	    if( cmp->compare(val,best_score) ) {
	      best_score = val;
	      tbk.init(mIdx,this,type,val);
	    }
	  }
	}
      }
    }
  }
  return pair<AnnotationScore,MatrixItemContents>(best_score,tbk);
}

static const char* 
extractCodonFrom3to5(const AnnotationItem& item, const DnaStr& theSequence, int bp_need, const MatrixIdx& mIdx) { 
  char* codon = justDontAskAboutThisOptimizationHackShit; 
  codon[0] = codon[1] = codon[2] = '\0';
  if(bp_need <= 0) return codon;
  //assert(side==dsu::eLeft);
  const MatrixItemContents& cell = item.getBest(mIdx);
  if(!cell.getPred())
    return codon;

  int start = cell.getPred()->getEnd5();
  int stop = item.getEnd3();
  if( mIdx.getSide() == dsu::eLeft) {
    const AnnotationItem* predCell2 = cell.getPred()->getBest(cell).getPred();
    if(predCell2 ) {
      start = predCell2->getEnd5();
    }
    stop = cell.getPred()->getEnd3();
  }
  int clen = stop-start+1, cnt = 0;

  while(bp_need > 0 && clen > 0) {
    --bp_need;
    codon[bp_need] = theSequence[(stop-cnt)];
    --clen;
    ++cnt;
  }
  if(bp_need == 0) return codon;
  //cout<<"Warning this is an exon of length: "<<start<<" "<<stop<<endl;
  return NULL;
}

bool 
AnnotationItem::checkOrf(const AnnotationItem* item, const MatrixIdx& mi, int end5, int end3, const DnaStr& str) {
  const int frame = mi.getFrame();
  int codon_len = frame>0?3-frame:frame;
  if( item && frame > 0 ) {
    char* codon = const_cast<char*>(extractCodonFrom3to5(*item,str,codon_len,mi));
    if( codon ) {
      if( frame==1 ) {
	codon[2] = str[end5];
      } else if(frame==2) {
	codon[1] = str[end5];
	codon[2] = str[end5+1];
      }
      if(str.isStop(codon))
	return false;
    }
  }
  int start = end5+frame;
  int stop = end3;
  int diff = ((stop-start)+1) % 3;
  stop-=diff;
  bool
#if _RUNTIME_DEBUG < 2
    val = noStops(str,start,stop);
#else
    val = noStops(str,start,stop,DEBUGON);
#endif
    return val;
}

void
AnnotationItem::update(const MatrixIdx& mi, const MatrixIdx& pIdx, int type, const AnnotationScore& is,
		       const AnnotationItem& litem) {
  const Options& opt = *Options::getInstance();
  const Compare* cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    cmp = GtCompare::getInstance();
  MatrixItemContents& mc = getBest(mi);
  if( AnnotationScore::badVal() != is && cmp->compare(is,mc.getScore()) ) {
    mc.init(pIdx,&litem,type,is);
  }
}

void
AnnotationItem::update(const MatrixIdx& mi, const MatrixIdx& pIdx, int type, const AnnotationScore& is) {
  const Options& opt = *Options::getInstance();
  const Compare* cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    cmp = GtCompare::getInstance();
  MatrixItemContents& mc = getBest(mi);
  if( AnnotationScore::badVal() != is && cmp->compare(is,mc.getScore()) ) {
    mc.init(pIdx,NULL,type,is);
  }
}

pair<int,dsu::Strand_t>
AnnotationItem::getPrevExonType(const MatrixIdx& mIdx) const {
  const Options& opt = *Options::getInstance();
  const Compare* cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    cmp = GtCompare::getInstance();
  const MatrixItemContents& pIdx = getBest(mIdx);
  int type = pIdx.getType();
  if(pIdx.getFrame()==-1)
    type = -1;
  dsu::Strand_t strnd = pIdx.getStrnd();
  return pair<int,dsu::Strand_t>(type,strnd);
}

void
AnnotationItem::dpUpdate(const AnnotationItem& litem, std::vector<scr_pr>& inv_scores, bool isContig, const CustomScore& cs) {
  SingleExon::update(isContig,litem,*this,inv_scores,cs);
  InitialExon::update(isContig,litem,*this,inv_scores,cs);
  InternalExon::update(isContig,litem,*this,inv_scores,cs);
  TerminalExon::update(isContig,litem,*this,inv_scores,cs);
  IntronInterval::update(litem,*this,inv_scores,cs);
  IgrInterval::update(litem,*this,inv_scores,cs);
}

void
  AnnotationItem::dpUpdate(std::vector<scr_pr>& inv_scores, bool isContig, bool allowPartialGene, const CustomScore& cs) {
  // IgrInterval and IntronInterval MUST be called BEFORE the exon intervals to set the "left hand side"
  IgrInterval::update_base(*this,inv_scores);

  if( allowPartialGene ) {
    IntronInterval::update_base(*this,inv_scores);
  }
  TerminalExon::update_base(isContig,*this,inv_scores,allowPartialGene,cs);
  InitialExon::update_base(isContig,*this,inv_scores,allowPartialGene,cs);
  SingleExon::update_base(isContig,*this,inv_scores,cs);
  if( allowPartialGene ) {
    InternalExon::update_base(isContig,*this,inv_scores,cs);
  }
}

int AnnotationItem::getEnd5AdjForStop(int model, dsu::Strand_t strnd) const {
  if( (model == Terminal::getInstance().getScoreIdx() || model == Single::getInstance().getScoreIdx()) && strnd == dsu::eNeg) {
    if( _negSeq.isStop(getEnd5()) ) return getEnd5();
    else {
      if( getEnd5()-3 >= 0 ) {
	if( _negSeq.isStop(getEnd5()-3 ) )
	  return getEnd5()-3;
      }
    }
    return getEnd5();
  }
  return getEnd5();
}

int AnnotationItem::getEnd3AdjForStop(int model, dsu::Strand_t strnd) const {
  if( (model == Terminal::getInstance().getScoreIdx() 
       || model == Single::getInstance().getScoreIdx()) && strnd == dsu::ePos) {
    if( _posSeq.isStop(getEnd3()-2) ) return getEnd3();
    else {
      if((getEnd3()+2 < (signed)_posSeq.length())) {
	if( _posSeq.isStop(getEnd3()+1 ) )
	  return getEnd3()+3;
      }
    }
    return getEnd3();
  }
  return getEnd3();
}

void
AnnotationItem::print(std::ostream& os, const CustomScore& cs) const {
  const AnnotationItem& ai = *this;
  os<<"-------------------------------------"<<endl;
  for(unsigned strnd = 0; strnd < 2; ++strnd) {
    for(unsigned frame = 0; frame < 3; ++frame) {
      for(unsigned side = 0; side < 2; ++side) {
	for(unsigned bnd = 0; bnd < 2; ++bnd) {
	  MatrixIdx mIdx(frame,(dsu::Strand_t)strnd,(dsu::side_t)side,(dsu::gbnd_t)bnd);
	  if(ai.getBest(mIdx).getType() != -1) {
	    os<<ai.getBest(mIdx)<<"|"<<mIdx<<endl;
	  }
	}
      }
    }
  }
  ai._seq.print(os,cs);
  os<<"-------------------------------------"<<endl;
}
