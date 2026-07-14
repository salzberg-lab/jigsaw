#include "SingleExon.h"
#include "MatrixIdx.h"
#include "genemodels.h"
#include "AnnotationItem.h"
#include "TransProb.h"
#include "SingleExon.h"
#include "CustomScore.h"
#include "ScoreInterval.h"

using namespace dsu;

// Single
void
SingleExon::update(bool isContig, const AnnotationItem& litem, AnnotationItem& item, std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  if(isContig) {
    update_pos(litem,item,inv_scores,cs);
  }
}
	     
void
SingleExon::update_pos(const AnnotationItem& litem, AnnotationItem& ritem, std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const TransProb& tp = TransProb::getInstance();
  const RegionPrediction& left = litem.getRegionPrediction();
  const RegionPrediction& rght = ritem.getRegionPrediction();
  const int cSNGL = Single::getInstance().getScoreIdx();
  dsu::Strand_t lstrnd=dsu::ePos,rstrnd=dsu::ePos;
  const int lframe = 0, rframe=0;
  const dsu::side_t lside=dsu::eLeft;
  const dsu::side_t rside=dsu::eRght;
  int length = (ritem.getEnd3()-litem.getEnd5())+1;
  // must be a multiple of 3!
  if( (length % 3) != 0 ) {
    return;
  }
  for(unsigned cand = 0; cand < 2; ++cand) {
    bool valdSngl=false;
    if(cand==0) {
      lstrnd=dsu::ePos,rstrnd=dsu::ePos;
      valdSngl = left.hasStart(lstrnd) && rght.hasStop(rstrnd);
    } else if(cand==1) {
      lstrnd=dsu::eNeg,rstrnd=dsu::eNeg;
      valdSngl = left.hasStop(lstrnd) && rght.hasStart(rstrnd);
    }
    if(!valdSngl) {
      continue;
    }
    const MatrixIdx rIdx(rframe,rstrnd,rside,eEndOfGene);
    const MatrixIdx lIdx(lframe,lstrnd,lside,eEndOfGene);
    if( litem.getBest(lIdx).getScore() != AnnotationScore::badVal()) {
      const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
      bool validEx = true;
      if(ptype.first != -1)
	validEx = tp.getProb(ptype.first,ptype.second,cSNGL,rstrnd)>0; 
      if( validEx ) {
	int loff=0,roff=3;
	if(lstrnd==dsu::eNeg) {
	  loff=3;
	  roff=0;
	}
	if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5()+loff,ritem.getEnd3()-roff,litem.getSeq(rstrnd)) ) {
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cSNGL,inv_scores[rstrnd].first,cs);
	  ritem.update(rIdx,lIdx,cSNGL,is,litem);
	}
      }
    }
  }
}

void
SingleExon::update_base(bool, AnnotationItem& item, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const int cSINGLE = Single::getInstance().getScoreIdx();
  const RegionPrediction& reg = item.getRegionPrediction();
  dsu::Strand_t strnd=dsu::ePos;
  if( item.length() % 3 != 0 ) {
    return;
  }
  if(reg.hasStart(strnd) && reg.hasStop(strnd)) {
    const int lframe = 0, rframe=0;
    const dsu::side_t lside=eLeft;
    const MatrixIdx lIdx(lframe,strnd,lside,eEndOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,reg.getEnd5(),reg.getEnd3()-3,item.getSeq(dsu::ePos))) {
      const dsu::side_t rside=eRght;
      const MatrixIdx rIdx(rframe,strnd,rside,eEndOfGene);
      const AnnotationScore is2=ScoreInterval::score(item,lIdx,item,rIdx,cSINGLE,0,cs);
      item.update(rIdx,lIdx,cSINGLE,is2,item);
    }
  } 
  strnd=dsu::eNeg;
  if(reg.hasStop(strnd) && reg.hasStart(strnd)) {
    const int lframe = 0, rframe=0;
    const MatrixIdx lIdx(lframe,dsu::eNeg,dsu::eLeft,eEndOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,reg.getEnd5()+3,reg.getEnd3(),item.getSeq(dsu::eNeg)) ) {
      const MatrixIdx rIdx(rframe,dsu::eNeg,dsu::eRght,eEndOfGene);
      const AnnotationScore is2=ScoreInterval::score(item,lIdx,item,rIdx,cSINGLE,0,cs);
      item.update(rIdx,lIdx,cSINGLE,is2,item);
    }
  } 
}
