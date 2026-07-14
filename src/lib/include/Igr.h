//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef IGR_H
#define IGR_H

#include "GeneModelType.h"

/**
 * <b>Responsibilities:</b><br>
 * instance of a <i>GeneModelType</i>.  Defines the behavior<br>
 * of a "single" exon.<br>
 * <br>
 * <b>Collaborator:</b><br>
 * GeneModelType 
 */
class Igr : public GeneModelType {
protected:

  Igr();

public:
  void print(std::ostream&,const string&) const;

  const int getScoreIdx() const { return 5; }
  const string& getTypeStr() const;
  static const Igr&  getInstance() { return _igr; };

private:
  static Igr _igr;

  /**
   * @link
   * @shapeType PatternLink
   * @pattern Igrton
   * @supplierRole Igrton factory 
   */
  /*# Igr _single1; */
};

#endif // IGR_H
