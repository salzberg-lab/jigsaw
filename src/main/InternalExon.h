#ifndef _INTERNALEXON_H
#define _INTERNALEXON_H

#include <vector>
#include "localtypedef.h"

class AnnotationItem;
class CustomScore;

class InternalExon {
 public:
  static void
    update(bool,const AnnotationItem&, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);

  static void
    update_base(bool, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore& );
};

#endif // _INTERNALEXON_H
