#ifndef _LOCALTYPEDEF_H
#define _LOCALTYPEDEF_H
#include "AnnotationScore.h"
#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 3 )
 #include <utility>
 using namespace std;
#else
 #include <pair.h>
#endif

//typedef pair<AnnotationScore,AnnotationScore> scr_pr; 
typedef pair<AnnotationScore, // protein coding interval
		       pair<AnnotationScore, // non-coding interval
				      AnnotationScore // intron interval score
					  > 
				> scr_pr; 

#endif //_LOCALTYPEDEF_Y
