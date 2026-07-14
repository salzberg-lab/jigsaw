//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef TERMINAL_H
#define TERMINAL_H

#include "GeneModelType.h"

/**
 * <b>Responsibilities:</b><br>
 * Instance of <i>GeneModelType</i>, defines the behavior for<br>
 * a "terminal" exon type.<br>
 * <br>
 * <b>Collaborators:</b><br>
 * GeneModelType 
 */
class Terminal : public GeneModelType {
protected:

  Terminal();

public:
  const int getScoreIdx() const { return 2; }
  const string& getTypeStr() const;
  void print(std::ostream&, const string&) const;
  static const Terminal& getInstance() { return _terminal; }

private:
  static Terminal _terminal;

  /**
   * @link
   * @shapeType PatternLink
   * @pattern Singleton
   * @supplierRole Singleton factory 
   */
  /*# Terminal _terminal1; */
};

#endif //TERMINAL_H
