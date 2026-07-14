#ifndef _TERMINALEXON_H
#define _TERMINALEXON_H

#include <vector>
#include "localtypedef.h"

class AnnotationItem;
class CustomScore;

class TerminalExon {
 public:
  static void
    update(bool,const AnnotationItem&, AnnotationItem&, const std::vector<scr_pr>&, const CustomScore&);

  static void
    update_base(bool, AnnotationItem&, const std::vector<scr_pr>&,bool, const CustomScore&);

};

#endif // _TERMINALEXON_H
