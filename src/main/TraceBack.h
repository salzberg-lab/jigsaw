#ifndef _TRACEBACK_H
#define _TRACEBACK_H

#include <list>
#include "AnnotationItem.h"
#include "GenePrediction.h"
#include "MatrixItemContents.h"

class TraceBack : public list<GenePrediction*> {
 public:
  //typedef AnnotationItem::strnd_frm_t strnd_frm_t;
  //typedef AnnotationItem::tback tback;

 public:
  TraceBack() { }

  void run(const MatrixItemContents&);

  pair<GenePrediction*,const MatrixItemContents*> 
  extractExon(const MatrixItemContents&,bool);
};

#endif //_TRACEBACK_H
