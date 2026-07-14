#ifndef _IGRINTERVAL_H
#define _IGRINTERVAL_H

#include <vector>
#include "localtypedef.h"

class AnnotationItem;
class CustomScore;

class IgrInterval {
 public:
  static void
    update(const AnnotationItem&, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);

  static void
    update_base(AnnotationItem&, const std::vector<scr_pr>&);

};

#endif // _IGRINTERVAL_H
