#ifndef _INITIALEXON_H
#define _INITIALEXON_H

#include <vector>
#include "localtypedef.h"

class AnnotationItem;
class CustomScore;

class InitialExon {
 public:
  static void
    update(bool,const AnnotationItem&, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);

  static void
    update_base(bool, AnnotationItem&, const std::vector<scr_pr>&,bool, const CustomScore&);

  static void
    update_neg(const AnnotationItem& , AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);

  static void
    update_pos(const AnnotationItem&, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);
};

#endif // _INITIALEXON_H
