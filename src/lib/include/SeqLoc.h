//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef SEQLOC_H
#define SEQLOC_H

#include <string>
#include <iostream>
#include <cstdlib>
#include <dsu.h>
#include "DnaStr.h"

/**
 * <b>Responsibilities:</b><br>
 * This is a general interface for dealing with specific<BR>
 * locations on a double stranded sequence<br>
 * SeqLoc is implemented specifically with the idea of 
 * traversing through a double stranded sequence. This means
 * that 5' < 3' always, and that the sequence location 
 * may refer to exclusively to one strand or both
 * depending on the constructor 
 * <b>Collaborators:</b> <br>
 * None. SeqLoc is the Base Class of a hierarchy<br>
 */
class SeqLoc {
  friend std::ostream& operator<<(std::ostream&,const SeqLoc&);
  inline friend bool operator<(const SeqLoc&, const SeqLoc&);
public:
  /**
   * a sequence location is made up of a strand, (positive or negative, or both)
   * and the 5' 3' location. 
   * 5' is from the positive strand perspective always, so 5' < 3' always. 
   * The 5' and 3' coordinates are inclusive
   * so that the length of the sequence is 3'-5'+1
   *
   * @param strand - dsu::ePos,dsu::eNeg or dsu::eEither
   * @param int - 5' location
   * @param int - 3' location
   */
   SeqLoc(dsu::Strand_t, int, int);

   /**
    * copy constructor
    */
   SeqLoc(const SeqLoc&);

   /**
    * object does not manage dynamic memory allocation
    */
   virtual ~SeqLoc() { }

   /**
    * return which strand this location refers to
    */
   inline dsu::Strand_t getStrnd() const;

   /**
    * get 5' end (always from positive perspective)
    */
   inline int getEnd5() const;

   /**
    * get 3' end (always from positive perspective)
    */
   inline int getEnd3() const;

   /**
    * @param int - refers to an idx to a base pair 
    * in the sequence
    * <br>
    * returns true if the idx falls within the 5',3'
    * coordinates of this sequence location
    */
   bool inRange(int) const;

   /**
    * sequence location length: 3'-5'+1
    */
   unsigned length() const;

   /**
    * if 5' == -1 returns false, true otherwise
    */
   bool isEnd5Set() const;

   /**
    * if 3' == -1 returns false, true otherwise
    */
   bool isEnd3Set() const;
   
   /**
    * returns true if both sequence locations
    * are on the same strand, and have the same
    * 5' and 3' values
    */
   bool operator==(const SeqLoc& sl) const {
     return _strnd == sl._strnd && _end5 == sl._end5 && sl._end3 == _end3;
   }
	
   /**
    * object identifier, useful, for polymorophic applications of SeqLoc
    */
   virtual const string& getTypeStr() const;

   /**
    * print sequence location information 
    * @param const string& - prints always take an optional string
    * argument that allows the user to prepend a string to the output
    */
   virtual void print(std::ostream&,const string& = "") const;

   /**
    * returns 5' (derived classes overide this behavior)
    */
   virtual int getClosestEnd5(int) const { return getEnd5(); } 

   /**
    * returns 3' (derived classes overide this behavior)
    */
   virtual int getClosestEnd3(int,int) const { return getEnd3(); } 

   /**
    * Precondition, int,int are already within the boundaries
    * of this sequence location this only applies to seqlocations
    * broken into sub sequence elements, e.g. DataSetAbstractProduct
    */
   virtual bool isActualOverlap(int,int) const { return true; }

   /**
    * user can set 5' end
    */
   void setEnd5(int);

   /**
    * user can set 3' end
    */
   void setEnd3(int);

   /**
    * user can set 3' end
    */
   void setStrnd(dsu::Strand_t strnd) { _strnd = strnd; }

   friend bool isOverlap(const SeqLoc&, const SeqLoc&);

private:
   dsu::Strand_t _strnd;
   int _end5, _end3;
};

#include "inline/SeqLoc.hpp"

#endif //SEQLOC_H
