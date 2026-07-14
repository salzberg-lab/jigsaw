#include "InitialExon.h"
#include "MatrixIdx.h"
#include "AnnotationItem.h"
#include "genemodels.h"
#include "TransProb.h"
#include "CustomScore.h"
#include "ScoreInterval.h"

using namespace dsu;

void
InitialExon::update(bool isContig, const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  if(isContig) {
    update_pos(litem,ritem,inv_scores,cs);
    update_neg(litem,ritem,inv_scores,cs);
  }
}

void
InitialExon::update_pos(const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const int INIT = Initial::getInstance().getScoreIdx();
  const TransProb& tp = TransProb::getInstance();
  const RegionPrediction& left = litem.getRegionPrediction();
  const RegionPrediction& rght = ritem.getRegionPrediction();
  const dsu::Strand_t lstrnd=dsu::ePos,rstrnd=dsu::ePos;
  const DnaStr& str = litem.getSeq(dsu::ePos);
  const int len = (rght.getEnd3()-left.getEnd5()) + 1;
  if( (len >= 3 && left.hasStart(lstrnd)) && (rght.isDonrBndry(rstrnd) || rght.getEnd3()==(signed)str.length()) ) {
    const int lframe = 0,rframe=(ritem.getEnd3()-litem.getEnd5()+1)%3;
    const dsu::side_t lside=dsu::eLeft;
    const dsu::side_t rside=dsu::eRght;
    const MatrixIdx lIdx(lframe,lstrnd,lside,eEndOfGene);
    const MatrixIdx rIdx(rframe,rstrnd,rside,eMdlOfGene);
    if( litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
      if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5(),ritem.getEnd3(),litem.getSeq(dsu::ePos)) ) {
	const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
	bool validExon = true;
	if(ptype.first != -1)
	  validExon = tp.getProb(ptype.first,ptype.second,INIT,dsu::ePos)>0;
	if(validExon) {
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,INIT,inv_scores[rstrnd].first,cs);
	  ritem.update(rIdx,lIdx,INIT,is,litem);
	}
      }
    }
  } 
}

void
InitialExon::update_neg(const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const TransProb& tp = TransProb::getInstance();
  const RegionPrediction& left = litem.getRegionPrediction();
  const RegionPrediction& rght = ritem.getRegionPrediction();
  const dsu::Strand_t lstrnd=dsu::eNeg,rstrnd=dsu::eNeg;
  const DnaStr& str = litem.getSeq(dsu::ePos);
  const int len = (rght.getEnd3()-left.getEnd5()) + 1;
  if((left.isDonrBndry(lstrnd) || left.getEnd5()==1) && (len >= 3 && rght.hasStart(rstrnd) || rght.getEnd3()==(signed)str.length())  ) {
    const dsu::side_t lside=dsu::eLeft;
    const dsu::side_t rside=dsu::eRght;
    const int lframe = (rght.getEnd3()-left.getEnd5()+1) % 3, rframe=0 ;
    const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
    const MatrixIdx rIdx(rframe,rstrnd,rside,eEndOfGene);
    if( litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
      if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5(),ritem.getEnd3(),ritem.getSeq(dsu::eNeg)) ) {
	const int INIT = Initial::getInstance().getScoreIdx();
	const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
	bool validExon = true;
	if(ptype.first != -1)
	  validExon = tp.getProb(ptype.first,ptype.second,INIT,dsu::eNeg)>0;
	if(validExon) {
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,INIT,inv_scores[rstrnd].first,cs);
	  ritem.update(rIdx,lIdx,INIT,is,litem);
	}
      }
    }
  } 
}

void
InitialExon::update_base(bool, AnnotationItem& item, const std::vector<scr_pr>& inv_scores, bool allowPartialGenes, const CustomScore& cs) {
  const int cINITIAL = Initial::getInstance().getScoreIdx();
  const RegionPrediction& reg = item.getRegionPrediction();
  const MatrixIdx offIdx(-1,dsu::eEither,dsu::eSide,eEndOfGene);
  dsu::Strand_t strnd = dsu::ePos;
  const int len = item.length();
  if(len >= 3 && reg.hasStart(strnd) && reg.isDonrBndry(strnd) && strnd==dsu::ePos) {
    const int lframe = 0;
    const int rframe = len % 3;
    const dsu::Strand_t lstrnd = dsu::ePos;
    const dsu::side_t lside=eLeft;
    const MatrixIdx lIdx(lframe,lstrnd,lside,eEndOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,reg.getEnd5(),reg.getEnd3(),item.getSeq(dsu::ePos)) ) {
      const dsu::side_t rside=eRght;
      const MatrixIdx rIdx(rframe,strnd,rside,eMdlOfGene);
      const AnnotationScore is2 = ScoreInterval::score(item,lIdx,item,rIdx,cINITIAL,0,cs);
      item.update(rIdx,lIdx,cINITIAL,is2,item);
    }
  } 
  strnd = dsu::eNeg;
  // Initial
  if(len >= 3 && reg.getEnd5()==1 && reg.hasStart(strnd) && strnd==dsu::eNeg && allowPartialGenes ) {
    const int lframe = reg.getEnd3() % 3, rframe=0 ;
    const MatrixIdx lIdx(lframe,dsu::eNeg,dsu::eLeft,eEndOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,item.getEnd5(),item.getEnd3(),item.getSeq(dsu::eNeg)) ) {
      const AnnotationScore is=inv_scores[strnd].first;
      const MatrixIdx rIdx(rframe,dsu::eNeg,dsu::eRght,eEndOfGene);
      const AnnotationScore is2=ScoreInterval::scorePartial1ToRight(item,rIdx,cINITIAL,is,cs);
      item.update(rIdx,offIdx,cINITIAL,is2);
    }
  } 
  strnd = dsu::eNeg;
  // Initial
  if(len >= 3 && reg.hasStart(strnd) && reg.isDonrBndry(strnd) && strnd==dsu::eNeg && allowPartialGenes) {
    const int rframe = 0;
    const int len = reg.getEnd3()-reg.getEnd5()+1;
    const int lframe = len % 3;
    const dsu::Strand_t lstrnd = dsu::eNeg;
    const dsu::side_t lside=eLeft;
    const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,item.getEnd5(),item.getEnd3(),item.getSeq(dsu::eNeg)) ) {
      const dsu::side_t rside=eRght;
      const MatrixIdx rIdx(rframe,strnd,rside,eEndOfGene);
      const AnnotationScore is2 = ScoreInterval::score(item,lIdx,item,rIdx,cINITIAL,0,cs);
      item.update(rIdx,lIdx,cINITIAL,is2,item);
    }
  } 
}
