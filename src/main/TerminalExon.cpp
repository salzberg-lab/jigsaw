#include "TerminalExon.h"
#include "MatrixIdx.h"
#include "AnnotationItem.h"
#include "genemodels.h"
#include "TransProb.h"
#include "ScoreInterval.h"

using namespace dsu;

void
TerminalExon::update(bool isContig, const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores,
		     const CustomScore& cs) {
  if(isContig) {
    const TransProb& tp = TransProb::getInstance();
    const int cTERM = Terminal::getInstance().getScoreIdx();
    const RegionPrediction& left = litem.getRegionPrediction();
    const RegionPrediction& rght = ritem.getRegionPrediction();
    const DnaStr& str = litem.getSeq(dsu::ePos);
    // terminal
    dsu::Strand_t strnd = dsu::ePos;
    bool leftIsValid = (left.isAccpBndry(strnd) || left.getEnd5()==1);
    bool rghtIsValid = (rght.hasStop(strnd) || rght.getEnd3()==(signed)str.length());
    if(isContig && leftIsValid && rghtIsValid) {
      const int lframe = (rght.getEnd3()-left.getEnd5()+1) % 3, rframe=0 ;
      const dsu::side_t lside=eLeft;
      const dsu::side_t rside=eRght;
      const MatrixIdx lIdx(lframe,strnd,lside,eMdlOfGene);
      const MatrixIdx rIdx(rframe,strnd,rside,eEndOfGene);
      const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
      bool validEx = true;
      if(ptype.first != -1)
	validEx = tp.getProb(ptype.first,ptype.second,cTERM,strnd) > 0;
      if( validEx && litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
	if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5(),ritem.getEnd3()-3,ritem.getSeq(dsu::ePos)) ) {
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cTERM,inv_scores[strnd].first,cs);
	  ritem.update(rIdx,lIdx,cTERM,is,litem);
	}
      }
    } 
    // terminal
    strnd = dsu::eNeg;
    leftIsValid = (left.hasStop(strnd) || left.getEnd5()==1);
    rghtIsValid = (rght.isAccpBndry(strnd) || rght.getEnd3()==(signed)str.length());
    if(isContig && leftIsValid && rghtIsValid) {
      dsu::side_t lside=eLeft;
      dsu::side_t rside=eRght;
      const int lframe = 0,rframe=(ritem.getEnd3()-litem.getEnd5()+1)%3;
      const MatrixIdx lIdx(lframe,strnd,lside,eEndOfGene);
      const MatrixIdx rIdx(rframe,strnd,rside,eMdlOfGene);
      const pair<int,dsu::Strand_t> ptype = litem.getPrevExonType(lIdx);
      bool validEx = true;
      if(ptype.first != -1)
	validEx = tp.getProb(ptype.first,ptype.second,cTERM,strnd) > 0;
      if( validEx && litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
	if( AnnotationItem::checkOrf(&litem,lIdx,litem.getEnd5()+3,ritem.getEnd3(),ritem.getSeq(dsu::eNeg)) ) {
	  const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cTERM,inv_scores[strnd].first,cs);
	  ritem.update(rIdx,lIdx,cTERM,is,litem);
	}
      }
    } 
  }
}

void
TerminalExon::update_base(bool isContig, AnnotationItem& item, const std::vector<scr_pr>& inv_scores, bool allowPartialGenes,
			  const CustomScore& cs) {
  // Terminal
  const int cTERMINAL = Terminal::getInstance().getScoreIdx();
  dsu::Strand_t strnd = dsu::eNeg;
  const RegionPrediction& reg = item.getRegionPrediction();
  if(reg.hasStop(strnd) && reg.isAccpBndry(strnd)) {
    const int lframe = 0;
    const int len = item.getEnd3()-item.getEnd5()+1;
    const int rframe = len % 3;
    const MatrixIdx lIdx(lframe,dsu::eNeg,dsu::eLeft,eEndOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,item.getEnd5()+3,item.getEnd3(),item.getSeq(dsu::eNeg)) ) {
      const MatrixIdx rIdx(rframe,dsu::eNeg,dsu::eRght,eMdlOfGene);
      const AnnotationScore is2=ScoreInterval::score(item,lIdx,item,rIdx,cTERMINAL,0,cs);
      item.update(rIdx,lIdx,cTERMINAL,is2,item);
    }
  } 
  // Terminal
  strnd = dsu::ePos;
  if(reg.isAccpBndry(strnd) && reg.hasStop(strnd) && allowPartialGenes) {
    const int len = item.getEnd3()-item.getEnd5()+1;
    const int lframe = len % 3, rframe=0 ;
    const dsu::Strand_t lstrnd = dsu::ePos;
    const dsu::side_t lside=eLeft;
    const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,item.getEnd5(),item.getEnd3()-3,item.getSeq(dsu::ePos)) ) {
      const dsu::side_t rside=eRght;
      const MatrixIdx rIdx(rframe,strnd,rside,eEndOfGene);
      const AnnotationScore is2=ScoreInterval::score(item,lIdx,item,rIdx,cTERMINAL,0,cs);
      item.update(rIdx,lIdx,cTERMINAL,is2,item);
    }
  } 
  strnd = dsu::ePos;
  // Terminal
  if(reg.hasStop(strnd) && reg.getEnd5()==1 && allowPartialGenes ) {
    const int lframe = item.getEnd3() % 3, rframe=0 ;
    const dsu::Strand_t lstrnd = dsu::ePos;
    const dsu::side_t lside=eLeft;
    const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
    if( AnnotationItem::checkOrf(NULL,lIdx,1,item.getEnd3()-3,item.getSeq(dsu::ePos)) ) {
      const AnnotationScore is=inv_scores[strnd].first;
      const dsu::side_t rside=eRght;
      const MatrixIdx rIdx(rframe,strnd,rside,eEndOfGene);
      const AnnotationScore is2=ScoreInterval::scorePartial1ToRight(item,rIdx,cTERMINAL,is,cs);
      item.update(rIdx,lIdx,cTERMINAL,is2);
    }
  }     
}
