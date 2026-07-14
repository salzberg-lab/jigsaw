//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
inline unsigned SeqLoc::length() const { return (abs(_end3-_end5)+1); }
inline dsu::Strand_t SeqLoc::getStrnd() const { return _strnd; }
inline int SeqLoc::getEnd5() const { return _end5; }
inline int SeqLoc::getEnd3() const { return _end3; }

bool operator<(const SeqLoc& as, const SeqLoc& bs) {
  if(as.getEnd5() < bs.getEnd5()) return true;
  if(as.getEnd5() == bs.getEnd5() &&
     as.getEnd3() < bs.getEnd3() ) return true;
  return false;
}
