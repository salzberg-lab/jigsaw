//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
bool AnnotationScore::operator<=(const AnnotationScore& score) const {
  return _val <= score._val;
}

bool AnnotationScore::operator>=(const AnnotationScore& score) const {
  return _val >= score._val;
}

bool AnnotationScore::operator>(const AnnotationScore& score) const {
  return _val > score._val;
}

bool AnnotationScore::operator==(const AnnotationScore& score) const {
  return _val == score._val;
}

AnnotationScore operator*(const AnnotationScore& score, double val) {
  return AnnotationScore( (score._val*val ));
}
