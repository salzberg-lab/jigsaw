//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef INTRON_H
#define INTRON_H

#include "GeneModelType.h"

/**
 * <b>Responsibilities:</b><br>
 * instance of a <i>GeneModelType</i>.  Defines the behavior<br>
 * of a "single" exon.<br>
 * <br>
 * <b>Collaborator:</b><br>
 * GeneModelType 
 */
class Intron : public GeneModelType {
protected:

  Intron();

public:
  void print(std::ostream&,const string&) const;

  const int getScoreIdx() const { return 4; }
  const string& getTypeStr() const;
  static const Intron&  getInstance() { return _intron; };

private:
  static Intron _intron;

  /**
   * @link
   * @shapeType PatternLink
   * @pattern Intronton
   * @supplierRole Intronton factory 
   */
  /*# Intron _single1; */
};

#endif // INTRON_H
