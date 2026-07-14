//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "SeqLoc.h"
#include <iostream>
#include <assert.h>
#include <algorithm>

SeqLoc::SeqLoc(dsu::Strand_t strnd, int end5, int end3) : _strnd(strnd), _end5(end5), _end3(end3) {
	if(_end5 > _end3) {
	   std::swap(_end5,_end3);
   }	
}

SeqLoc::SeqLoc(const SeqLoc& sl) : _strnd(sl._strnd), _end5(sl._end5), _end3(sl._end3) {
}

const string internalStr = "SeqLoc";
const string& SeqLoc::getTypeStr() const {
  return internalStr;
}

bool SeqLoc::inRange(int idx) const {
  return idx >= getEnd5() && idx <= getEnd3();
}

void SeqLoc::setEnd5(int val) { 
  _end5 = val;
}

void SeqLoc::setEnd3(int val) { 
  _end3 = val;
}

bool SeqLoc::isEnd3Set() const {
  return _end3==-1 ? false : true;
}

bool SeqLoc::isEnd5Set() const {
  return _end5==-1 ? false : true;
}

void SeqLoc::print(std::ostream& os, const string& str2) const {
  os<<str2<<strnd2str(getStrnd())<<" "<<getEnd5()<<" "<<getEnd3(); 
} 

std::ostream& operator<<(std::ostream& os, const SeqLoc& sl) {
  sl.print(os);
  return os;
}

bool
isOverlap(const SeqLoc& aloc, const SeqLoc& bloc) {
  return ( aloc.inRange(bloc.getEnd5()) || aloc.inRange(bloc.getEnd3()) ||
	   bloc.inRange(aloc.getEnd5()) || bloc.inRange(aloc.getEnd3())); 
}
