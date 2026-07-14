#include "InternalExon.h"
#include "MatrixIdx.h"
#include "AnnotationItem.h"
#include "genemodels.h"
#include "TransProb.h"
#include "ScoreInterval.h"
#include "CustomScore.h"

using namespace dsu;
#if _RUNTIME_DEBUG >= 3
extern bool DEBUGON;
#endif

void
InternalExon::update(bool isContig, const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  if(!isContig) {
    return;
  }
  const TransProb& tp = TransProb::getInstance();
  const RegionPrediction& left = litem.getRegionPrediction();
  const RegionPrediction& rght = ritem.getRegionPrediction();
  const int cINTERNAL = Internal::getInstance().getScoreIdx();
  const DnaStr& str = litem.getSeq(dsu::ePos);
  // internal
  for(unsigned lframe = 0; lframe < 3; ++lframe) {
    dsu::Strand_t lstrnd=dsu::ePos,rstrnd=dsu::ePos;
    const dsu::side_t lside=dsu::eLeft;
    const dsu::side_t rside=dsu::eRght;
    const int len = ritem.getEnd3()-litem.getEnd5()+1;
    const int rframe = abs((int)(len-lframe)) % 3;
    bool leftIsValid = (left.isAccpBndry(lstrnd) || left.getEnd5()==1);
    bool rghtIsValid = (rght.isDonrBndry(rstrnd) || rght.getEnd3()==(signed)str.length());
    if(leftIsValid && rghtIsValid) {
      const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
      const MatrixIdx rIdx(rframe,rstrnd,rside,eMdlOfGene);
      const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
      const bool validEx = tp.getProb(ptype.first,ptype.second,cINTERNAL,rstrnd) > 0;
      if( validEx && litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
	if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5(),ritem.getEnd3(),litem.getSeq(dsu::ePos))) {
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cINTERNAL,inv_scores[rstrnd].first,cs);
	  ritem.update(rIdx,lIdx,cINTERNAL,is,litem);
	}
      }
    }
    lstrnd=dsu::eNeg;
    rstrnd=dsu::eNeg;
    leftIsValid = (left.isDonrBndry(lstrnd) || left.getEnd5()==1);
    rghtIsValid = (rght.isAccpBndry(rstrnd) || rght.getEnd3()==(signed)str.length());
    if(leftIsValid && rghtIsValid) {
      const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
      const MatrixIdx rIdx(rframe,rstrnd,rside,eMdlOfGene);
      const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
      const bool validEx = tp.getProb(ptype.first,ptype.second,cINTERNAL,rstrnd) > 0;
      if( validEx && litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
	if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5(),ritem.getEnd3(),litem.getSeq(dsu::eNeg)) ) {
#if _RUNTIME_DEBUG >= 3
	  if( DEBUGON ) {
	    int stop = 0;
	  }
#endif
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cINTERNAL,inv_scores[rstrnd].first,cs);
	  ritem.update(rIdx,lIdx,cINTERNAL,is,litem);
	}
      }
    }
  } 
}

void
InternalExon::update_base(bool isContig, AnnotationItem& item, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const int cINTERNAL = Internal::getInstance().getScoreIdx();
  dsu::Strand_t strnd = dsu::ePos;
  const RegionPrediction& reg = item.getRegionPrediction();
  // this handles running off the left of the sequence
  // internal
  if(reg.getEnd5() == 1 && reg.isDonrBndry(strnd)) {
    for(unsigned lframe = 0; lframe < 3; ++lframe) {
      const int len = reg.getEnd3();
      const int rframe = abs((int)(len-lframe)) % 3;
      const dsu::side_t lside=eLeft;
      const MatrixIdx lIdx(lframe,strnd,lside,eMdlOfGene);
      if( AnnotationItem::checkOrf(NULL,lIdx,1,reg.getEnd3(),item.getSeq(dsu::ePos)) ) {
	const AnnotationScore is=inv_scores[strnd].first;
	const MatrixIdx rIdx(rframe,dsu::eNeg,dsu::eRght,eMdlOfGene);
	const AnnotationScore is2=ScoreInterval::scorePartial1ToRight(item,rIdx,cINTERNAL,is,cs);
	item.update(rIdx,lIdx,cINTERNAL,is2);
      }
    }
  }     
  // internal
  if(reg.isDonrBndry(strnd) && reg.isAccpBndry(strnd) ) {
    const int len = reg.getEnd3()-reg.getEnd5()+1;
    for(unsigned lframe = 0; lframe < 3; ++lframe) {
      const int rframe = abs((int)(len-lframe)) % 3;
      const dsu::Strand_t lstrnd = dsu::ePos;
      const dsu::side_t lside=eLeft;
      const dsu::side_t rside=eRght;
      const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
      if( AnnotationItem::checkOrf(NULL,lIdx,reg.getEnd5(),reg.getEnd3(),item.getSeq(dsu::ePos)) ) {
	const MatrixIdx rIdx(rframe,strnd,rside,eMdlOfGene);
	const AnnotationScore is2=ScoreInterval::score(item,lIdx,item,rIdx,cINTERNAL,0,cs);
	item.update(rIdx,lIdx,cINTERNAL,is2,item);
      }
    }
  } 
  strnd = dsu::eNeg;
  // Internal
  if(reg.getEnd5()==1 && reg.isAccpBndry(strnd) ) {
    for(unsigned lframe = 0; lframe < 3; ++lframe) {
      const int len = reg.getEnd3();
      const int rframe = abs((int)(len-lframe)) % 3;
      const MatrixIdx lIdx(lframe,dsu::eNeg,dsu::eLeft,eMdlOfGene);
      if( AnnotationItem::checkOrf(NULL,lIdx,1,reg.getEnd3(),item.getSeq(dsu::eNeg)) ) {
	const AnnotationScore is=inv_scores[strnd].first;
	const MatrixIdx rIdx(rframe,dsu::eNeg,dsu::eRght,eMdlOfGene);
	const AnnotationScore is2=ScoreInterval::scorePartial1ToRight(item,rIdx,cINTERNAL,is,cs);
	item.update(rIdx,lIdx,cINTERNAL,is2);
      }
    }
  } 
  // Internal
  if(reg.isDonrBndry(strnd) && reg.isAccpBndry(strnd) ) {
    const int len = reg.getEnd3()-reg.getEnd5()+1;
    for(unsigned lframe = 0; lframe < 3; ++lframe) {
      const int rframe = abs((int)(len-lframe)) % 3;
      const MatrixIdx lIdx(lframe,dsu::eNeg,dsu::eLeft,eMdlOfGene);
      if( AnnotationItem::checkOrf(NULL,lIdx,reg.getEnd5(),reg.getEnd3(),item.getSeq(dsu::eNeg)) ) {
	const MatrixIdx rIdx(rframe,strnd,dsu::eRght,eMdlOfGene);
	const AnnotationScore is2=ScoreInterval::score(item,lIdx,item,rIdx,cINTERNAL,0,cs);
	item.update(rIdx,lIdx,cINTERNAL,is2,item);
      }
    }
  }
}
