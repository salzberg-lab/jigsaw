//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEMODELS_H
#define GENEMODELS_H

#include "Initial.h"
#include "Terminal.h"
#include "Single.h"
#include "Internal.h"
#include "Igr.h"
#include "Intron.h"
#include "dsu.h"

namespace genemodels {

  const string& idx2str(int);
  /**
   ** input "right side" gene model type, 
   ** "left side" gene model type, bool is 
   ** contiguous sequence
  **/

  bool isNewExonInModel(int, dsu::Strand_t, int, dsu::Strand_t, bool );
  bool isInSameGeneModel(int, dsu::Strand_t, int, dsu::Strand_t, bool );
  inline bool isContiguous(int end5, int end3);
  inline bool isContiguous(int end5, dsu::Strand_t, int , dsu::Strand_t);

}

#include "inline/genemodels.hpp"

#endif
