//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef _EXONDISTR_H
#define _EXONDISTR_H

#include <vector>
#include <string>
using std::string;

/**
 * @class ExonDistr
 * 
 * singleton, records exon prob. values
 * The probability of N exons in a gene
 * 
 * if a "null" file name is passed
 * then default is all gene model
 * exon number distributions are
 * equally likely
 */

class ExonDistr {
 public:
  double getProb(int idx) const;
 public:
  static void createInstance(const string&);
  static const ExonDistr& getInstance();
 protected:
  ExonDistr(const string&);

 private:
  std::vector<double> _pdistr;
  /**
   * the purpose of the minimum factor
   * is compute a reasonable estimate for
   * exon probabilities that have
   * not been directly observed, rather
   * than just setting them to 0.
   * Instead, find out what the minimum
   * probability value is, and the
   * order of magnitude it is at (10^-x)
   * where x = _minFactor
   * then 1/exonnumber * minFactor
   * is the estimate, this approach 
   * assumes that more exons per gene model
   * are less likely
   */
  double _minFactor;
 private:
  static ExonDistr* _instance;
};

#endif //_EXONDISTR_H
