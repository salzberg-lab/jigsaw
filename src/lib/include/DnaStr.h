//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef DNASTR_H
#define DNASTR_H

#include <string>
#include <dsu.h>
#include <strstream>
#include <stdexcept>
#include <assert.h>

/**
 * Facilitates some useful operations on a particular double stranded dna sequence 
 * (e.g. positive or negative)
 */

class DnaStr {
 public:
  /**
   * simply maintains a const reference to the double strand sequence 
   * and sets an internal variable to remember which strand to conduct
   * operations on.
   */
  DnaStr(const string& str, dsu::Strand_t strnd) : _str(str), _strnd(strnd) { }

  /**
   * get const reference to double strand
   */
  const string& getStr() const { return _str; }

  /**
   * gc %
   * @param - int 5' end (default = 1)
   * @param - int 3' end (default = sequence length)
   */
  double gcPcnt(int = 1, int = -1) const;

  /**
   * query which side of the strand the object is set to operate on 
   */
  dsu::Strand_t getStrnd() const { return _strnd; }

  /**
   * access a single base at a specified location in the string
   */
  const char operator[](int idx) const
#if _RUNTIME_DEBUG < 3
    {
    assert(idx >= 1 && (idx-1) < (signed)_str.length());
    return _strnd==dsu::eNeg ? dsu::comp(_str[ idx-1]) : _str[idx-1];
  }
#else
  ;
#endif

  /**
   * query length of strand
   */
  inline unsigned length() const { return _str.length(); }

  /**
   * check for a start codon
   * sorry, this should be made static
   * @param const char[3] = input codon
   */
  inline bool isStart(const char [3]) const;

  /**
   * check for a start codon begining at a specified location
   * in the strand 
   */
  inline bool isStart(int)const;

  /**
   * check for a stop codon
   * sorry, this should be made static
   * @param const char[3] = input codon
   */
  inline bool isStop(const char [3]) const;

  /**
   * check for stop codon starting a specified locationi in strand
   */
  inline bool isStop(int bp) const; 

  /**
   * print 5' to 3'
   */
  friend std::ostream& operator<<(std::ostream& os, const DnaStr& str) {
		for(unsigned i = 0; i < str._str.length(); ++i)
			os<<str[i];
		return os;
  }

#if _RUNTIME_DEBUG >= 1
  void printCodon(int start) const {
    if(_strnd == dsu::ePos)
      std::cout<<operator[](start)<<operator[](start+1)<<operator[](start+2)<<std::endl;
    else
      std::cout<<operator[](start+2)<<operator[](start+1)<<operator[](start)<<std::endl;
  }
#endif
private:
  inline static bool isU(char);
  inline static bool isA(char);
  inline static bool isG(char);
  inline static bool isC(char);
  inline static bool isR(char);

 private:
  const string& _str;
  dsu::Strand_t _strnd;

};

#include "inline/DnaStr.hpp"

#endif // DNASTR_H
