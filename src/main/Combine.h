//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef _COMBINE_H
#define _COMBINE_H

#include "Annotation.h"
#include "InputData.h"

class RegionPrediction;
class CustomScore;

/**
 * <b>Purpose:</b><br>
 * Initiates the main combiner algorithm.  This is the public interface to the combiner algor.<br> 
 */
class Combine {
 public:
	/**<b>Input:</b> <br>
	 * Cronstructing this object, runs the combiner algorithm
	 * @param InputData& : reference to all the gene analysis information
	 * which has been organized along the dna sequence by the InputData object<br>
	 * @param const string& : a const reference to the actual sequence<br>
	 * to look.  Algorithm operates on both strands at the same time.<br>
	 * @param const CustomScore& = reference to the user defined weighting function
	 * (see the public api for a defintion)
	 * @param int - bp cutoff, when to stop comparing combinations, when x number of bp's sperate
	 * the two coding regions.
	 * <b>PostCond:</b> <br>
	 * An a scoring of all posible combinations, and a recording of the best combination
	 * (which must be traced back through)
	 **/
  Combine(const InputData &, const string&, bool,const CustomScore&, int);

  ~Combine();

  void run();

	
  list<RegionPrediction*>& getOutput() { return _lst; }
  const list<RegionPrediction*>& getOutput() const { return _lst; }
  bool empty() const;
  Annotation& getAnnotation() { return *_annotation; }

 private:
  void getAllRelaventData(int,InputData::const_iterator&,bool);
  void mergeScoreInsert();
  void score(list<RegionPrediction*>&);

 private:

    /**
     * @supplierCardinality 1 
     */
  const InputData& _inData;
  const string& _theSequence;
  Annotation* _annotation;
  const CustomScore* _cstmScr;
  list<RegionPrediction*> _lst;
};

#endif // _COMBINE_H
