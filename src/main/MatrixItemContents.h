#ifndef _MATRIXITEMCONTENTS_H
#define _MATRIXITEMCONTENTS_H

#include "MatrixIdx.h"
//#include "AnnotationItem.h"
#include "AnnotationScore.h"

class AnnotationItem;

class MatrixItemContents : public MatrixIdx {
 public:
  MatrixItemContents() : _pred(NULL), _type(-1) { }
  MatrixItemContents(const MatrixIdx& mi, const AnnotationItem* pred, int type, const AnnotationScore& scr) : 
    MatrixIdx(mi), _pred(pred), _type(type), _score(scr)   
    {}
  MatrixItemContents(const MatrixItemContents& mi) : MatrixIdx(mi), _pred(mi._pred), _type(mi._type), _score(mi._score) 
    {}

  void init(const MatrixIdx& mi, const AnnotationItem* pred, int type, const AnnotationScore& scr);
  void init(const MatrixItemContents&);

  AnnotationScore getScore() const { return _score;}
  const AnnotationItem* getPred() const { return _pred; }
  int getType() const { return _type; }

  friend std::ostream&
  operator<<(std::ostream&,const MatrixItemContents&);

 private:
  const AnnotationItem* _pred;
  int _type;
  AnnotationScore _score;
};


#endif // _MATRIXITEMCONTENTS_H
