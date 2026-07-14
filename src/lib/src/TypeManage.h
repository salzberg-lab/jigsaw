#ifndef _TYPEMANAGE_H
#define _TYPEMANAGE_H

#include "MatchPrediction.h"
#include "GenePrediction.h"
#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 3 )
 #include <utility>
 using namespace std;
#else
 #include <pair.h>
#endif

class TypeManage {
 public:
  static const MatchPrediction* 
    checkForMatchPrediction(const Prediction*, int); 

  static pair<const Prediction*,const Prediction*>
    checkForIntronPrediction(const Prediction*, int, const Prediction*, int); 

  static pair<const Prediction*,const Prediction*>
    checkForIntronPrediction(const Prediction*, int, int); 

  static const MatchPrediction* 
    checkForMatchPrediction(const Prediction*, int, int); 

  static const MatchPrediction* 
    checkForMatchPrediction5(const Prediction*, int); 

  static const MatchPrediction* 
    checkForMatchPrediction3(const Prediction*, int); 

  static const GenePrediction* 
    checkForGenePrediction(const Prediction*, int); 

  static const GenePrediction* 
    checkForGenePrediction(const Prediction*, int, int); 

  static const GenePrediction* 
    checkForGenePrediction5(const Prediction*, int); 

  static const GenePrediction* 
    checkForGenePrediction3(const Prediction*, int); 
};
#endif // _TYPEMANAGE_H
