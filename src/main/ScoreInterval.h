#ifndef _SCOREINTERVAL_H
#define _SCOREINTERVAL_H

#include "AnnotationScore.h"
#include "dsu.h"

class AnnotationItem;
class RegionPrediction;
class MatrixIdx;
class CustomScore;

class ScoreInterval {
 public:
  static AnnotationScore 
    score(const AnnotationItem&, const MatrixIdx&, 
	  const AnnotationItem&, const MatrixIdx&, 
	  const int, const AnnotationScore&, const CustomScore&);
  
  static AnnotationScore 
    scorePartial1ToRight(const AnnotationItem& rItem, const MatrixIdx& rIdx,
			 const int type, const AnnotationScore&, const CustomScore&);

  static AnnotationScore 
    scorePartialToEnd(const AnnotationItem&, const MatrixIdx&,
		      const int, const AnnotationScore&); 

  static AnnotationScore 
    intronPenalty(const int,  const int, const AnnotationScore&); 

  static AnnotationScore 
    scoreOne(const AnnotationItem&, const AnnotationItem&, const int, dsu::Strand_t, const AnnotationScore&, const CustomScore&);

 private:
  static AnnotationScore 
    igrPenalty(const RegionPrediction&, const RegionPrediction&, const AnnotationScore&); 

  static AnnotationScore 
    intronPenalty(const RegionPrediction&, const RegionPrediction&, const AnnotationScore&); 

};

#endif // _SCOREINTERVAL_H
