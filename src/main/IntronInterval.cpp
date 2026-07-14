#include "IntronInterval.h"
#include "ParamFile.h"
#include "ScoreInterval.h"
#include "CustomScore.h"
#include "AnnotationItem.h"
#include "genemodels.h"
#include "TransProb.h"
#include "MatrixIdx.h"
#include "ScoreInterval.h"

using std::cout;
using std::endl;
using namespace dsu;

#if _RUNTIME_DEBUG >= 2
extern bool DEBUGON;
#endif


static char justDontAskAboutThisOptimizationHackShit[3];
static const char* 
extractCodonFrom3to5(const AnnotationItem& item, const DnaStr& theSequence, int bp_need, const MatrixIdx& mi) {
  char* codon = justDontAskAboutThisOptimizationHackShit; 
  codon[0] = codon[1] = codon[2] = '\0';
  if(bp_need <= 0) return codon;
  const AnnotationItem* predCell = item.getBest(mi).getPred();
  if(!predCell)
    return codon;

  const MatrixIdx& pIdx = item.getBest(mi).getIdx();
  int start = predCell->getEnd5();
  int stop = item.getEnd3();
  if( mi.getSide() == dsu::eLeft) {
    const AnnotationItem* predCell2 = predCell->getBest(pIdx).getPred();
    if(predCell2 ) {
      start = predCell2->getEnd5();
    }
    stop = predCell->getEnd3();
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


void
IntronInterval::update(const AnnotationItem& litem, AnnotationItem& ritem, const std::vector<scr_pr>& inv_scores, const CustomScore& cs) {
  const int cINTRON = Intron::getInstance().getScoreIdx();
  const TransProb& tp = TransProb::getInstance();
  const RegionPrediction& left = litem.getRegionPrediction();
  const RegionPrediction& rght = ritem.getRegionPrediction();
  const ParamFile& pf = *(ParamFile::getInstance());
  const int minIntronLen=pf.getMinIntronLength();
  // intron
  const dsu::side_t lside = dsu::eRght;
  const dsu::side_t rside = dsu::eLeft;
  for(unsigned rframe = 0; rframe < 3; ++rframe) {
    dsu::Strand_t lstrnd = dsu::ePos, rstrnd = dsu::ePos;
    const int lframe = rframe > 0 ? (3-(rframe % 3)) : rframe;
    if(left.isDonrBndry(lstrnd) && rght.isAccpBndry(rstrnd) && ((rght.getEnd5()-left.getEnd3())>minIntronLen)) {
      const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
      const MatrixIdx rIdx(rframe,rstrnd,rside,eMdlOfGene);
      const int codon_len = 3;
      const DnaStr& str = litem.getSeq(lstrnd);
      bool noStop=true;
      if(lframe == 0) {
	char* codon = const_cast<char*>(extractCodonFrom3to5(litem,str,codon_len,lIdx));
	if( codon ) {
	  noStop = !str.isStop(codon);
	}
      }
      if(noStop) {
	const int ltype = litem.getBest(lIdx).getType();
	const bool validExon = tp.getProb(ltype,lstrnd,cINTRON,rstrnd)>0;
	if(validExon) {
	  if( litem.getBest(lIdx).getScore() != AnnotationScore::badVal() && noStop ) {
	    const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cINTRON,inv_scores[rstrnd].second.second,cs);
	    ritem.update(rIdx,lIdx,cINTRON,is,litem);
	  }
	}
      }
    }
    lstrnd = dsu::eNeg;
    rstrnd = dsu::eNeg;
    // intron
    if(left.isAccpBndry(lstrnd) && rght.isDonrBndry(rstrnd) && ((rght.getEnd5()-left.getEnd3())> minIntronLen) ) {
      const MatrixIdx lIdx(lframe,lstrnd,lside,eMdlOfGene);
      const MatrixIdx rIdx(rframe,rstrnd,rside,eMdlOfGene);
      const DnaStr& str = litem.getSeq(lstrnd);
      bool noStop=true;
      if( rght.length() >=3 && rframe==0 ) {
	noStop = !str.isStop(ritem.getEnd5());
      }
      if(noStop) {
	const int ltype = litem.getBest(lIdx).getType();
	const bool validExon = tp.getProb(ltype,lstrnd,cINTRON,rstrnd)>0;
	if(validExon) {
	  if( litem.getBest(lIdx).getScore() != AnnotationScore::badVal() ) {
	    const AnnotationScore is=ScoreInterval::score(litem,lIdx,ritem,rIdx,cINTRON,inv_scores[rstrnd].second.second,cs);
	    ritem.update(rIdx,lIdx,cINTRON,is,litem);
	  }
	}
      }
    }
  }
}

void
IntronInterval::update_base(AnnotationItem& item, const std::vector<scr_pr>& inv_scores) {
  const RegionPrediction& rght = item.getRegionPrediction();
  const int cINTRON = Intron::getInstance().getScoreIdx();
  const MatrixIdx offIdx(-1,dsu::eEither,dsu::eSide,eEndOfGene);
#if _RUNTIME_DEBUG >= 2
  if( item.getEnd5() == 639) {
    int stop=0;
  }
#endif
  for(unsigned  i = 0; i < 2; ++i) {
    const dsu::Strand_t rstrnd = (dsu::Strand_t)i;
    /**
     * Here's the deal, I'm not sure the best way to handle partial
     * genes when it comes to intron evidence right now gene prediction
     * programs are the only evidence source which are capable of predicting
     * a partial intron (splice predictions too but not implimented).
     * for this reason I want to use the IGR score, since in
     * the representations, an intron will never be predicted and the
     * way to keep this fair for now is to use the same score
     * here the intron score is commented out
    const AnnotationScore is=inv_scores[rstrnd].second.second;
    **/
    //const AnnotationScore is=ScoreInterval::intronPenalty(0,item.getEnd5(),inv_scores[rstrnd].second.first);
    const AnnotationScore is=ScoreInterval::intronPenalty(0,item.getEnd5(),inv_scores[rstrnd].second.second);
    // Intron
    if(rght.isAccpBndry(rstrnd) && rstrnd==dsu::ePos) {
      for(unsigned rframe = 0; rframe < 3; ++rframe) {
	const dsu::side_t rside=eLeft;
	const MatrixIdx lIdx(rframe,rstrnd,rside,eMdlOfGene);
	item.update(lIdx,offIdx,cINTRON,is);
      }
    } 
    // Intron
    if(rght.isDonrBndry(rstrnd) &&rstrnd==dsu::eNeg ) {
      for(unsigned lframe = 0; lframe < 3; ++lframe) {
	const MatrixIdx lIdx(lframe,dsu::eNeg,dsu::eLeft,eMdlOfGene);
	item.update(lIdx,offIdx,cINTRON,is);
      }
    } 
  }
}
