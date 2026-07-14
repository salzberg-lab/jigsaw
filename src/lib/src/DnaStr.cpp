//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "DnaStr.h"

#if _RUNTIME_DEBUG >= 3
const char DnaStr::operator[](int idx) const {
  assert(idx >= 1 && (idx-1) < (signed)_str.length());
  return _strnd==dsu::eNeg ? dsu::comp(_str[ idx-1]) : _str[idx-1];
}
#endif

double DnaStr::gcPcnt(int end5, int end3) const {
  end3 = end3==-1?length():end3;
  int cnt = 0;
  for(int iter = end5; iter < end3; ++iter) {
    char base = operator[](iter);
    if( isG(base) || isC(base) ) {
      ++cnt;
    }
  }
  return (double)cnt/(double)(end3-end5+1);
}
