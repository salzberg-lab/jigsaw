//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEMODELTYPE_H
#define GENEMODELTYPE_H

#include <string>
using std::string;

/**
 * <b>Purpose:</b><br>
 * Interface - defines the behavior of an idividual exon
 * model type. e.g Terminal, Initial, etc.<br>
 * Important note that GeneModelTypes are really not used to improve gene prediction at the moment
 * at the moment.  Currently elements are considered
 * coding, and model infromation is not used, if at sometime
 * important information becomes available to incorporate gene model 
 * information, we'll be ready...
 * @stereotype interface 
 */
class GeneModelType {
 public:
  virtual void print(std::ostream&, const string&) const = 0;

  /**
   * <b>Postcond:</b><br>
   * returns a string identifier that identifies type (for printing)<br> 
   */
  virtual const string& getTypeStr() const = 0;

  /**
   * <b>Postcond:</b><br>
   * returns a unique integer id representing the model type that can be<br>
   * used to access a unique location in an array.<br> 
   */
  virtual const int getScoreIdx() const  = 0;
};

#endif // GENEMODELTYPE_H
