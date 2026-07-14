//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef ANNOTATIONSCORE_H
#define ANNOTATIONSCORE_H

#include <iostream>

//class ostream;
//class istream;

/**
 * <b>Purpose:</b><br>
 * Abstracts the notion of a "score", has ended up behaving<br>
 * just like a "double", but hey how was I suppose to know!<br> 
 */
class AnnotationScore {
 public:
  /**
   * Initializes to "AnnotationScore::badVal()"
   */
  AnnotationScore() : _val(AnnotationScore::badVal().getVal()) { }
  AnnotationScore(double val) : _val(val) {}
  AnnotationScore(int val) : _val(static_cast<double>(val)) {}
  AnnotationScore(const AnnotationScore& as) : _val(as._val) { }
 ~AnnotationScore() { }

/**
 * Checks to see if score is equivalent to badVal()
 */
 bool isValid() const { return badVal()._val != _val; }

/**
 * const value defaultScore() == 1.0
 */
 inline static AnnotationScore defaultScore() {
  return AnnotationScore(1.0);
 }

/**
 * const value = -1.0
 */
inline static AnnotationScore minVal() {
  return AnnotationScore(-1.0);
}

/**
 * const value = -1.0
 */
inline static AnnotationScore badVal() {
  return AnnotationScore(-1.0);
}
    
inline friend AnnotationScore operator*(const AnnotationScore&,double);
inline bool operator>(const AnnotationScore&) const;
inline bool operator<=(const AnnotationScore&) const;
inline bool operator>=(const AnnotationScore&) const;
inline bool operator==(const AnnotationScore&) const;
bool operator!=(const AnnotationScore& as) const {
  return (_val != as._val);
}
inline AnnotationScore& operator=(const AnnotationScore& in) {
  _val = in._val;
  return *this;
}
inline AnnotationScore& operator/=(const AnnotationScore& in) {
  _val /= in._val;
  return *this;
}
inline AnnotationScore& operator-=(const AnnotationScore& in) {
    _val -= in._val;
    return *this;
}
inline AnnotationScore& operator+=(const AnnotationScore& in) {
    _val += in._val;
    return *this;
}
inline AnnotationScore& operator*=(const AnnotationScore& in) {
    _val *= in._val;
    return *this;
}

 friend inline AnnotationScore operator*(const AnnotationScore& a, const AnnotationScore& b) {
   return a._val*b._val;
 } 
friend bool operator<(const AnnotationScore& ain, const AnnotationScore& bin) {
  const double aVal = ain.getVal();
  const double bVal = bin.getVal();
  const bool result = (aVal < bVal);
  return result;
}
inline friend AnnotationScore operator-(const AnnotationScore& as1, const AnnotationScore& as2) {
    const AnnotationScore rVal = AnnotationScore(as1._val-as2._val);
    return rVal;
  }

inline friend AnnotationScore operator+(const AnnotationScore& as1, const AnnotationScore& as2) 
    { return AnnotationScore( as1._val+as2._val ); }

inline friend AnnotationScore operator/(const AnnotationScore& as1, const AnnotationScore& as2) {
   return AnnotationScore(as1._val/as2._val); 
}

 inline friend AnnotationScore combineScore(const AnnotationScore& as1, const AnnotationScore& as2) {
    return (as1 + as2); 
} 


 friend std::istream& operator>>(std::istream&,AnnotationScore&);
 friend std::ostream& operator<<(std::ostream&,const AnnotationScore&);
 inline double getVal() const { return _val; }

private:
 double _val;

};

#include "inline/AnnotationScore.hpp"

#endif //ANNOTATIONSCORE_H
