//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef INITIAL_H
#define INITIAL_H

#include "GeneModelType.h"

/**
 * <b>Responsibilities:</b><br>
 * define behavior of an "initial" exon model.<br>
 * instance of the <i>GeneModelType</i><br>
 * <br>
 * <b>Collaborators:</b><br>
 * GeneModelType 
 */
class Initial : public GeneModelType {
protected:

  Initial();

public:
  void print(std::ostream& os, const string&) const;
  const string& getTypeStr() const;
  const int getScoreIdx() const { return 0; }

  static const Initial&  getInstance() { return _initial; };

private:
  static Initial _initial;

  /**
   * @link
   * @shapeType PatternLink
   * @pattern Singleton
   * @supplierRole Singleton factory 
   */
  /*# Initial _initial1; */
};


#endif //INITIAL_H
