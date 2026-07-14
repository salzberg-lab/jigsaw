#ifndef _EVIDENCETYPEREGISTRY_H
#define _EVIDENCETYPEREGISTRY_H

#include "EvidenceType.h"
#include <string>


/**
 * class ormatRegistry
 * The purpose of this class is to
 * prevent any format specific details
 * from permeating the combiner application
 * code.  The details for various format
 * handlers are kept hidden the development 
 * library.
 */

class EvidenceTypeRegistry {
 public:
  /** 
   * @param const string& : string identifier specifying a prediction type: acc,don,start,stop,coding,intron
   * @param const string& : string identifier specifying a program type: geneprediction, homology, spliceprediction
   * @param const string& : string identifying source, usually the filename
   * @param int : int identifier to make source unique (look into streamlining this)
   */
  static EvidenceType* 
    createNewEvidenceType(const std::string&, const std::string&, const std::string&, int, const std::string&, const double = 1.0);

  static EvidenceType* 
    createNewEvidenceTypeLinearBackwardCompatible(const std::string&, const std::string&, double, const std::string&, int, const std::string&);
};

#endif //_EVIDENCETYPEREGISTRY_H
