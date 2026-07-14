//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "genemodels.h"

// just an optimization 
static const string initialStr = "initial";
static const string internalStr = "internal";
static const string terminalStr = "terminal";
static const string singleStr = "single";
static const string errorStr = "error";

const string& genemodels::idx2str(int idx) {
  if(idx == Initial::getInstance().getScoreIdx())
      return initialStr;
  else if(idx == Internal::getInstance().getScoreIdx())
      return internalStr;
  else if(idx == Terminal::getInstance().getScoreIdx())
    return terminalStr;
  else if(idx == Single::getInstance().getScoreIdx())
    return singleStr;
  else
      return errorStr;
}

bool genemodels::isInSameGeneModel(int typePost, dsu::Strand_t currStrnd, int typePred, dsu::Strand_t predStrnd, bool isContig) {
  if(typePred == typePost && isContig && currStrnd==predStrnd) return true;
  if(currStrnd != predStrnd) return false;
  if(currStrnd == dsu::ePos) {
    if(typePred == Terminal::getInstance().getScoreIdx() || typePred == Single::getInstance().getScoreIdx()) return false;
  } else {
    if(typePred == Initial::getInstance().getScoreIdx() || typePred == Single::getInstance().getScoreIdx()) return false;
  }
  return true;
}

bool genemodels::isNewExonInModel(int typePost, dsu::Strand_t currStrnd, int typePred, dsu::Strand_t predStrnd, bool isContig) {
  return (isInSameGeneModel(typePost,currStrnd,typePred,predStrnd,isContig) && !isContig);
}
