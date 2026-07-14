//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef EVSCORES_H
#define EVSCORES_H

#include "AnnotationScore.h"
#include "dsu.h"
#include <assert.h>
#include <vector>
#include <list>
#include <iostream>
class AnnotationItem;
class RegionPrediction;

class EvScores {
 public:
  /** 
   * Defines different scorable conditions
   */
  enum Bnd_t { eAcc = 0,  /**< acceptor and start exon boundaries */
	       eDon,      /**< donor and stop exon boundaries */
	       eCoding,  /** non boundary conditions */
	       eStart,
	       eStop,
	       eIntron,
	       eNumBnd ,
	       eNonCoding,
	       eAccH = 0x8,
	       eDonH = 0x10,
	       eCodingH = 0x20,
	       eStartH = 0x40,
	       eStopH = 0x80,
	       eIgr = 0x100,
	       eIntr = 0x200,
	       eErr = 0x400
  };
  friend std::ostream& operator<<(std::ostream&, const EvScores&);

 public:
  EvScores();
  EvScores(const EvScores&);
  EvScores& operator=(const EvScores&);
  const AnnotationScore& operator[](int idx) const {
    assert(idx >= 0 && idx < eNumBnd);
    return _score[idx];
  }
  AnnotationScore& operator[](int idx) {
    assert(idx >= 0 && idx < eNumBnd);
    return _score[idx];
  }

 private:
  AnnotationScore _score[eNumBnd];

};

#endif // EVIDENCESCORES_H

