//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef SINGLE_H
#define SINGLE_H

#include "GeneModelType.h"

/**
 * <b>Responsibilities:</b><br>
 * instance of a <i>GeneModelType</i>.  Defines the behavior<br>
 * of a "single" exon.<br>
 * <br>
 * <b>Collaborator:</b><br>
 * GeneModelType 
 */
class Single : public GeneModelType {
protected:

  Single();

public:
  void print(std::ostream&,const string&) const;

  const int getScoreIdx() const { return 3; }
  const string& getTypeStr() const;
  static const Single&  getInstance() { return _single; };

private:
  static Single _single;

  /**
   * @link
   * @shapeType PatternLink
   * @pattern Singleton
   * @supplierRole Singleton factory 
   */
  /*# Single _single1; */
};

#endif // SINGLE_H
