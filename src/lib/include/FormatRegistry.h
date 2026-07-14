#ifndef _FORMATREGISTRY_H
#define _FORMATREGISTRY_H

#include "DataSetAbstractFactory.h"
#include <iostream>


/**
 * class FormatRegistry
 * The purpose of this class is to
 * prevent any format specific details
 * from permeating the combiner application
 * code.  The details for various format
 * handlers are kept hidden the development 
 * library.
 */

class FormatRegistry {
 public:
  /** 
   * @param const string& : string identifier specifying a file format
   * @param const string& : string identifier specifying a prediction type
   * recognized file format identifiers:
   * default
   * gff
   * btab
   * glimmerm, fgenesh, phat, genscan
   * genemarkhmm
   *
   * recognized prediction types:
   * geneprediction
   * homology
   * spliceprediction
   */
  static DataSetAbstractFactory* getFactory(const string&, const string&);

  static void printFormats(std::ostream&);
};

#endif //_FORMATREGISTRY_H
