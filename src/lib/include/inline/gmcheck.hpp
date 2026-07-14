//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
/** speed things up a bit... **/
static int INTERNAL_INT_ID;
static int INITIAL_INT_ID;
static int TERMINAL_INT_ID;
static int SINGLE_INT_ID;

static bool checkForAdjacentStopCodon(const RegionPrediction& predReg, const RegionPrediction& itemReg, 
				      const DnaStr& theSequence) {
  if(itemReg.getEnd5()-predReg.getEnd3() > 3) {
    char codon[3];
    int idx = predReg.getEnd3()+1;
    codon[0] = theSequence[idx];
    codon[1] = theSequence[idx+1];
    codon[2] = theSequence[idx+2];
    return theSequence.isStop(codon);
  }
  return false;
}

static bool isValidFrameFromLen(int start, int stop, int bp_need, int frame) {
  const int raw_len = stop-start+1;
  if(raw_len == 1) {
    if(bp_need == 1 && frame == 0) return true;
    if(bp_need == 2 && frame == 1) return true;
    if(bp_need == 3 && frame == 2) return true;
  } else if(raw_len == 2) {
    if(bp_need == 1 && frame == 2) return true;
    if(bp_need == 2 && frame == 0) return true;
    if(bp_need == 3 && frame == 1) return true;
  } else {
    int leading_bp = (raw_len-bp_need) % 3;
    int predFrame = AnnotationItem::leadingBp2Frame(leading_bp);
    if(predFrame == frame) return true;
  }
  return false;
}

static char justDontAskAboutThisOptimizationHackShit[3];
static const char* extractCodonFrom3to5(const AnnotationItem& item, const UnitLoc& ul,
					const DnaStr& theSequence, int bp_need) {
  char* codon = justDontAskAboutThisOptimizationHackShit; 
  if(bp_need <= 0) return codon;
  const int start = item.getEnd5(), stop = item.getEnd3();
  int clen = stop-start+1, cnt = 0;

  if(!isValidFrameFromLen(start,stop,bp_need,ul.getFrame())) return NULL;
  while(bp_need > 0 && clen > 0) {
    --bp_need;
    codon[bp_need] = theSequence[(stop-cnt)];
    --clen;
    ++cnt;
  }

  if(bp_need == 0) return codon;

  assert(ul.getStrnd()==theSequence.getStrnd());
  const UnitScore& unit = item.getMatrix().getUnitScore(ul);
  const AnnotationItem* predCell = unit.getPredCell();
  if(predCell) {
     const AnnotationItem& predItem = *predCell;
     const UnitLoc& pul = unit.getPredUnitLoc();
     return extractCodonFrom3to5(predItem,pul,theSequence,bp_need);
  }
  return NULL;
}

/** 
 ** We need to know if 'item' is part of a codon that is starting off a new model
 ** this only for frames 0,1
*/
static bool checkFirstCodonOfNewModel(const AnnotationItem& item, int itemType, int currFrame,
                   const AnnotationItem& pred, const UnitLoc& pl) {
  int predType = pl.getType();
  int predFrame = pl.getFrame();
  bool rghtBnd = pl.isRghtBnd();
  bool isBeginParse = pl.isBeginParse();
  //preCond
  assert(currFrame==0||currFrame==1);
  assert((item.length()==1&&currFrame==1) || (currFrame==0&&(item.length()==1||item.length()==2)));
  if( itemType == Internal::getInstance().getScoreIdx() || predType == Internal::getInstance().getScoreIdx()) 
    return false;
  if( item.getStrnd() == dsu::ePos && (itemType==Terminal::getInstance().getScoreIdx())) return false;
  if( item.getStrnd() == dsu::eNeg && (itemType==Initial::getInstance().getScoreIdx())) return false;

  bool isContig = genemodels::isContiguous(pred.getEnd3(),pred.getStrnd(),item.getEnd5(),item.getStrnd());
  // just lookoing to see if this is part of the 1st codon of a new model

  // if this is first frame, then this must be the
  // the start of the codon and so the predecessor has to be part 
  // of another exon
  if( currFrame==0 && !isContig ) return true;

  // if this is part of the frame 1
  // then the predecessor must be of length 1 and it must be the start
  // of the new model
  bool predIsStart = false;
  int exonCnt = 0;
  const UnitScore& us = pred.getMatrix().getUnitScore(UnitLoc(predType,predFrame,exonCnt,pred.getStrnd(),rghtBnd,isBeginParse));
  if(us.getPredCell()) {
    const AnnotationItem& predPred = *us.getPredCell();
    predIsStart = !genemodels::isContiguous(predPred.getEnd3(), us.getPredStrnd(),pred.getEnd5(),pred.getStrnd());
  } else {
    predIsStart = true;
  }
  if( currFrame==1 && pred.length()==1 && predIsStart ) return true;

  // don't worry about yet
  // currFrame==3
  return false;
}

/**
 **
 */
static bool passSpecialLookAhead(const AnnotationItem& item, int itemType, int currFrame,
                   const AnnotationItem& pred, const UnitLoc& pul, bool newModel) {
  // PreCond
  int predType = pul.getType();
  assert((item.length()==1&&currFrame==1) || (currFrame==0&&(item.length()==1||item.length()==2)));
  bool isContig = genemodels::isContiguous(pred.getEnd3(),pred.getStrnd(),item.getEnd5(),item.getStrnd());
  const DnaStr& str = item.getSeq();

  // even if contig, if the contig is only 2bp, can assume that we can check next bp, cuase an exon has to at least be 3bp
  if( (isContig && pred.length()==1) || !isContig) {
    // make sure that pred's pred is not contig, otherwise it is a >2bp exon and can't assume next bp will be included
    const UnitScore& us = pred.getMatrix().getUnitScore(pul);
    if(us.getPredCell() || !isContig) {
      bool checkCodon = !isContig;
      if(us.getPredCell() && isContig) {
	const AnnotationItem& predPred = *us.getPredCell();
	checkCodon = !genemodels::isContiguous(predPred.getEnd3(),us.getPredStrnd(),pred.getEnd5(),pred.getStrnd());
      }
      char codon[3];
      bool isStop = false;
      if(checkCodon && currFrame == 1) {
      // check that we're not at end of sequence
	if(item.getEnd3()+1 >= (signed)item.getSeq().length()) return true;
	assert(us.getStrnd()==pred.getStrnd());
	codon[0] = str[pred.getEnd3()];
	codon[1] = str[item.getEnd5()];
	codon[2] = str[item.getEnd5()+1];
	isStop = str.isStop(codon);
      } else if(checkCodon && currFrame == 0) {
	// skip at the end
	if(item.getEnd5()+2 >= (signed)str.length()) return true;
	codon[0] = str[item.getEnd5()];
	codon[1] = str[item.getEnd5()+1];
	codon[2] = str[item.getEnd5()+2];
	isStop = str.isStop(codon);
      }
      if(checkCodon) {
	if( !newModel ) {
	  if(item.getStrnd()==dsu::ePos && 
	     // model could end hear in which case this would be a valid terminal or single ending on a stop
	     (predType != Terminal::getInstance().getScoreIdx() || predType != Single::getInstance().getScoreIdx()) && 
	     (itemType != Terminal::getInstance().getScoreIdx() || itemType != Single::getInstance().getScoreIdx()) && isStop) {
	    return false;
	  } 
	  else if(item.getStrnd()==dsu::eNeg && isStop)
	    return false;
	} else if( newModel ) {
	  // if neg strnd, can start off with stop, so nothing to check
	  assert(itemType != Internal::getInstance().getScoreIdx());
	  if(dsu::ePos==item.getStrnd())
	    assert(itemType == Initial::getInstance().getScoreIdx() || itemType == Single::getInstance().getScoreIdx());
	  else if(dsu::eNeg==item.getStrnd()) {
	    assert(itemType == Terminal::getInstance().getScoreIdx() || itemType == Single::getInstance().getScoreIdx());
	  }
	  // positives start off with an an initial so we can check for stop
	  bool isAdjStop = checkForAdjacentStopCodon(pred.getRegionPrediction(),item.getRegionPrediction(),str);
	  if(item.getStrnd() == dsu::ePos && isStop) {
	    return false;
	    // negatives start off with an an terminal so we can check for no stop
	  } else if( item.getStrnd() == dsu::eNeg && (!isStop && !isAdjStop) )
	    return false;
	}
      }
    }
  }
  // made it through all the special checks
  return true;
}

static bool predEndsModel(const AnnotationItem& predItem, const UnitLoc& pul, int len) {
  const UnitScore& pc = predItem.getMatrix().getUnitScore(pul);
  if( len == predItem.length() ) return true;
  else if(pc.getPredCell()) {
    const AnnotationItem& prev = *pc.getPredCell();
    if( prev.getStrnd() != predItem.getStrnd() ||
       !genemodels::isContiguous(prev.getEnd3(),pc.getPredStrnd(),predItem.getEnd5(),predItem.getStrnd()) ) {
      return false;
    } else { 
      const UnitLoc& npc = pc.getPredUnitLoc();
      return predEndsModel(prev,npc,len-predItem.length());
    }
  } else {
    return false;
  }
}


/**
 ** The purpose of this procedure is to insure that the predecessor is a valid 
 ** interval that the current interval can be connected to.
 */
bool doProcess(const AnnotationItem& item, const AnnotationItem& predItem, const UnitLoc& pul) {
  bool isContig =  genemodels::isContiguous(predItem.getEnd3(), predItem.getStrnd(), item.getEnd5(), item.getStrnd());
  if(isContig) return true;

  // bp 3
  /* check if the were connecting with a predecessor that was scored/valid
  **/
  const UnitScore& unit = predItem.getMatrix().getUnitScore(pul);
  if(!unit.isValid()) return false;
  if(predItem.length() >= 3) return true;

  // bp 2
  const AnnotationItem* predPredIter = unit.getPredCell();
  if(!predPredIter) return false;
  const AnnotationItem& predPredItem = *predPredIter;
  int end3 = predPredItem.getEnd3();
  // contiguous
  if(end3+1 != predItem.getEnd5()) return false;
  if( predPredItem.length() + predItem.length() >= 3) return true;

  // bp 1
  const UnitScore& nextUnit = 
    predPredItem.getMatrix().getUnitScore(unit.getPredUnitLoc());
  if(!nextUnit.isValid()) return false;

  const AnnotationItem* nextPredIter = nextUnit.getPredCell();
  if(!nextPredIter) return false;
  const AnnotationItem& nextPredItem = *nextPredIter;
  int nextEnd3 = nextPredItem.getEnd3();
  if(nextEnd3+1 != end3) return false;
  return true;
}

static int currentFrame2NumNeed(int frm) {
  switch(frm) {
  case 2: return 2;
  case 1: return 1;
  case 0: return 0;
  default:
    assert(0);
    return -1;
  }
}

static int getProperFrameForCurrFrameOne(int len) {
  int predFrame = -1;
  const int curr_frame = 0;
  if(len>=3) {
    int leading_bp = (len-currentFrame2NumNeed(curr_frame)) % 3;
    predFrame = AnnotationItem::leadingBp2Frame(leading_bp);
  } else if(len == 2) {
    predFrame = 1;
  } else {
    predFrame = 2;
  }
  return predFrame;
}

#if _RUNTIME_DEBUG >= 1
static bool noStops(const DnaStr& seq, int start, int stop, bool verbose) {
#else
static bool noStops(const DnaStr& seq, int start, int stop) {
#endif 
  for(int i = start; i <= stop-2; i+=3) {
#if _RUNTIME_DEBUG >= 1
    if(verbose)
      seq.printCodon(i);
#endif
    if(seq.isStop(i)) {
      return false;
    }
  }
  return true;
}

static bool hasAdjStop(const AnnotationItem& curr, const AnnotationItem& pred) {
  assert(curr.getStrnd() == dsu::eNeg);
  if((curr.getEnd5()-pred.getEnd3()) >= 4) {
    const DnaStr& str = curr.getSeq();
    if(str.isStop(curr.getEnd5()-3)) return true;
  }
  return false;
}

static bool check4PredContig(const AnnotationItem& item, const UnitLoc& ul, int bp_need) {
  const int start = item.getEnd5(), stop = item.getEnd3();
  int clen = stop-start+1;
  bp_need -= clen;
  if( bp_need <= 0 ) return true;
  const UnitScore& unit = item.getMatrix().getUnitScore(ul);
  const AnnotationItem* predCell = unit.getPredCell();
  const UnitLoc& pul = unit.getPredUnitLoc();
  if(predCell) {
     const AnnotationItem& predItem = *predCell;
     bool isContig = genemodels::isContiguous(predItem.getEnd3(),pul.getStrnd(),item.getEnd5(),item.getStrnd());
     if( !isContig && bp_need == 1) {
       return false;
     } else
       return check4PredContig(predItem,pul,bp_need);
  }
  return false;
}

inline static bool ignoreHelp(int itemPost, int itemPre, dsu::Strand_t strnd, bool isContig) {
  // ignore gene model starting with "internal"
  if(itemPost == TERMINAL_INT_ID && itemPre == TERMINAL_INT_ID && !isContig) return true;
  if(itemPost == INITIAL_INT_ID && itemPre == INITIAL_INT_ID && !isContig) return true;

  if(strnd == dsu::ePos) {
    if(itemPost == SINGLE_INT_ID &&
       (itemPre == INITIAL_INT_ID || itemPre == INTERNAL_INT_ID)) return true;
    if(itemPost == TERMINAL_INT_ID &&
      itemPre == SINGLE_INT_ID) return true;
    if(itemPost == INTERNAL_INT_ID &&
     (itemPre != INTERNAL_INT_ID && itemPre != INITIAL_INT_ID))
    return true;
    if(itemPost == INITIAL_INT_ID && itemPre == INTERNAL_INT_ID) return true;
  } else if(strnd == dsu::eNeg) {
    if(itemPost == SINGLE_INT_ID &&
       (itemPre == TERMINAL_INT_ID || itemPre == INTERNAL_INT_ID)) return true;
    if(itemPost == INITIAL_INT_ID &&
      itemPre == SINGLE_INT_ID) return true;
    if(itemPost == INTERNAL_INT_ID &&
       (itemPre != INTERNAL_INT_ID && itemPre != TERMINAL_INT_ID))
    return true;
    if(itemPost == TERMINAL_INT_ID && itemPre == INTERNAL_INT_ID) return true;
  } else assert(0);
     
  return false;
}

inline static  bool ignoreConnection(int currType, dsu::Strand_t currStrnd, int predType, dsu::Strand_t predStrnd, bool isContig) {
  if(currStrnd == predStrnd && isContig) return false;

  if(currStrnd == dsu::ePos && predStrnd == dsu::ePos) {
    return ignoreHelp(currType,predType,currStrnd,isContig);
  } else if(currStrnd == dsu::ePos && predStrnd == dsu::eNeg) {
    if((currType==SINGLE_INT_ID||currType==INITIAL_INT_ID) && 
       (predType==SINGLE_INT_ID||predType==INITIAL_INT_ID)) return false;
    return true;
  } else if(currStrnd == dsu::eNeg && predStrnd == dsu::ePos) {
    if((currType==SINGLE_INT_ID||currType==TERMINAL_INT_ID) && 
       (predType==SINGLE_INT_ID||predType==TERMINAL_INT_ID)) return false;
    return true;
  } else if(currStrnd == dsu::eNeg && predStrnd == dsu::eNeg) {
    return ignoreHelp(currType,predType,currStrnd,isContig);
  }
  return true;
}

/**
 ** looking at the case where item is part of a start codon
 */
static bool isEndOfStartCodonPass(const AnnotationItem& item, int itemType, int currFrame,
				  const AnnotationItem& pred, int predType, int predFrame,
				  const DnaStr& theSequence) {
  return true;
}
