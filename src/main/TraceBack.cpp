#include "TraceBack.h"
#include "genemodels.h"
#include "AnnotationItem.h"

using std::cout;
using std::endl;

void
TraceBack::run(const MatrixItemContents& start) {
  MatrixItemContents ct(start);
  bool isFirst=true;
  while(ct.getPred()) {
    pair<GenePrediction*,const MatrixItemContents*> val = extractExon(ct,isFirst);
    isFirst=false;
    GenePrediction* exon = val.first;
    if(exon) {
      push_front(exon);
    }
    if(val.second) {
      ct = *val.second;
    } else {
      break;
    }
  }
}

/**
 * PreCondition: cell,unit points to a valid right side of an exon
 *               otherwise NULL is returned (i know not a PreCond..)
 **/
pair<GenePrediction*,const MatrixItemContents*> 
TraceBack::extractExon(const MatrixItemContents& cback, bool isFirst) {
  const AnnotationItem* rght = cback.getPred();
  int exon_type = cback.getType();
  const MatrixItemContents& tbk = rght->getBest(cback);
  const AnnotationItem* left = tbk.getPred();
#if _RUNTIME_DEBUG >= 2
  cout<<"Right: "<<rght->getEnd5()<<" "<<rght->getEnd3()<<endl;
  if( rght->getEnd5() == 747 && rght->getEnd3()== 870 ) {
    int stop = 0;
  }
#endif
  if(tbk.getType()== -1) {
    return pair<GenePrediction*,const MatrixItemContents*>(NULL,NULL);
  }
  if( !left ) {
    left = rght;
  }
  //only in the first cell of trace back does cback contain the exon type
  //in all other cases it should be contianed in the rght node should
  //is the right end of an exon, in the first case it's computed at the
  // end
  if( !isFirst ) {
    exon_type = tbk.getType();
  }
  int end5 = left->getEnd5AdjForStop(exon_type,tbk.getStrnd());
#if _RUNTIME_DEBUG >= 2
    cout<<"Left: "<<left->getEnd5()<<" "<<left->getEnd3()<<endl;
#endif
  int end3 = rght->getEnd3AdjForStop(exon_type,cback.getStrnd());
  const GeneModelType* exon_label = NULL;
  if( exon_type == Internal::getInstance().getScoreIdx() ) 
    exon_label = &(Internal::getInstance());
  else if( exon_type == Initial::getInstance().getScoreIdx() ) 
    exon_label = &(Initial::getInstance());
  else if( exon_type == Single::getInstance().getScoreIdx() ) 
    exon_label = &(Single::getInstance());
  else if( exon_type == Terminal::getInstance().getScoreIdx() ) 
    exon_label = &(Terminal::getInstance());
  assert(exon_label);
  int lframe = tbk.getFrame();
//#warning "Make sure this new allocation is properly delete"
  GenePrediction* exon = new GenePrediction(SeqLoc(cback.getStrnd(),end5,end3),1.0,-1,"combiner",*exon_label, lframe);
  const MatrixItemContents& tmp = left->getBest(tbk);
  return pair<GenePrediction*,const MatrixItemContents*>(exon,&tmp);
}
