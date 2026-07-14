#ifndef _INTRONINTERVAL_H
#define _INTRONINTERVAL_H

#include <vector>
#include "localtypedef.h"

class AnnotationItem;
class PrintInterval;
class RegionPrediction;
class CustomScore;

class IntronInterval {
 public:
  static void
    update(const AnnotationItem&, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);

  static void
    update_base(AnnotationItem&, const std::vector<scr_pr>&);

  //static void
    //train_print(const RegionPrediction&, const RegionPrediction&, const PrintInterval& );

};

#endif // _INTRONINTERVAL_H
