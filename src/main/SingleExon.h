#ifndef _SINGLEEXON_H
#define _SINGLEEXON_H

#include <vector>
#include "localtypedef.h"

class AnnotationItem;
class CustomScore;

class SingleExon {
 public:
  static void
    update(bool,const AnnotationItem&, AnnotationItem&, std::vector<scr_pr>&, const CustomScore&);
  static void
    update_base(bool isContig, AnnotationItem& item, const std::vector<scr_pr>& inv_scores, const CustomScore&);
  static void
    update_neg(const AnnotationItem& , AnnotationItem&, std::vector<scr_pr>&, const CustomScore&);
  static void
    update_pos(const AnnotationItem&, AnnotationItem&, std::vector<scr_pr>&, const CustomScore&);
};

#endif // _SINGLEEXON_H
