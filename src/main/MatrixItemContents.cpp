#include "MatrixItemContents.h"
#include "AnnotationItem.h"
#include "genemodels.h"

void
MatrixItemContents::init(const MatrixIdx& mi, const AnnotationItem* pred, int type, const AnnotationScore& scr) {
  MatrixIdx::init(mi);
  _pred = pred;
  _type = type;
  _score = scr;
}

static string
type2str(int typ) {
  string val;
  if( typ == Initial::getInstance().getScoreIdx()) {
    val = "init";
  } else if( typ == Single::getInstance().getScoreIdx()) {
    val = "sngl";
  } else if( typ == Internal::getInstance().getScoreIdx()) {
    val = "intr";
  } else if( typ == Terminal::getInstance().getScoreIdx()) {
    val = "term";
  } else if( typ == Intron::getInstance().getScoreIdx()) {
    val = "tron";
  } else if( typ == Igr::getInstance().getScoreIdx()) {
    val = "igr";
  }
  return val;
}

std::ostream&
operator<<(std::ostream& os, const MatrixItemContents& mi) {
  string invTypeStr = type2str(mi.getType());
  os<<(MatrixIdx)mi<<" "<<invTypeStr<<" "<<mi.getScore();
  if(mi.getPred()) {
    os<<" "<<mi.getPred()->getEnd5()<<" "<<mi.getPred()->getEnd3();
  }
  return os;
}
