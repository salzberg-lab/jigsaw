#include "IgrInterval.h"
#include "ParamFile.h"
#include "AnnotationItem.h"
#include "CustomScore.h"
#include "genemodels.h"
#include "TransProb.h"
#include "MatrixIdx.h"
#include "ScoreInterval.h"

using namespace dsu;

void
IgrInterval::update(const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const int cIGR = Igr::getInstance().getScoreIdx();
  const TransProb& tp = TransProb::getInstance();
  const RegionPrediction& left = litem.getRegionPrediction();
  const RegionPrediction& rght = ritem.getRegionPrediction();
  const ParamFile& pf = *(ParamFile::getInstance());
  const int minIntergenicLen=pf.getMinIgrLength();
  const bool igrLongEnough = (rght.getEnd5()-left.getEnd3())>minIntergenicLen;
  // intergenic
  if( (RegionPrediction::is_igr(left,rght,dsu::ePos,0) || RegionPrediction::is_igr(left,rght,dsu::eNeg,0)) && 
      igrLongEnough) {
    int lframe = 0, rframe = 0;
    const dsu::side_t lside = dsu::eRght;
    const dsu::side_t rside = dsu::eLeft;
    for(unsigned i = 0; i < 2; ++i) {
      const dsu::Strand_t lstrnd = (dsu::Strand_t)i;
      const MatrixIdx lIdx(lframe,lstrnd,lside,eEndOfGene);
      const int ltype = litem.getBest(lIdx).getType();
      if(ltype != -1) {
	for(unsigned j = 0; j < 2; ++j) {
	  const dsu::Strand_t rstrnd = (dsu::Strand_t)j;
	  const MatrixIdx rIdx(rframe,rstrnd,rside,eEndOfGene);
	  const bool validExon = tp.getProb(ltype,lstrnd,cIGR,rstrnd)>0;
	  if(validExon) {
	    if( litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
	      const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cIGR,inv_scores[rstrnd].second.first,cs);
	      ritem.update(rIdx,lIdx,cIGR,is,litem);
	    }  
	  }
	}
      } 
    }
  }
}

void
IgrInterval::update_base(AnnotationItem& item, const std::vector<scr_pr>& inv_scores) {
  const RegionPrediction& rght = item.getRegionPrediction();
  const int cIGR = Igr::getInstance().getScoreIdx();
  const MatrixIdx offIdx(-1,dsu::eEither,dsu::eSide,eEndOfGene);
  for(unsigned  i = 0; i < 2; ++i) {
    const dsu::Strand_t rstrnd = (dsu::Strand_t)i;
    // Intergenic
    if((rght.hasStart(rstrnd) && rstrnd==dsu::ePos) ||
       (rght.hasStop(rstrnd) && rstrnd==dsu::eNeg)) {
      const int lframe = 0;
      const AnnotationScore is=inv_scores[rstrnd].second.first;
      const dsu::side_t lside=eLeft;
      const MatrixIdx lIdx(lframe,rstrnd,lside,eEndOfGene);
      item.update(lIdx,offIdx,cIGR,is);
    } 
  }
}
