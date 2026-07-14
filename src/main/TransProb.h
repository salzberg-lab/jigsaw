#ifndef _TRANSPROB_H
#define _TRANSPROB_H

#include "dsu.h"
#include <vector>

class TransProb {
 public:
  TransProb();
  double
  getProb(int,dsu::Strand_t,int,dsu::Strand_t) const;

  static const TransProb& 
  getInstance();

  static void 
  create();

 private:
  void init();
  double& 
  setProb(int,dsu::Strand_t,int,dsu::Strand_t);
  
  std::vector<std::vector<double> > _vec;

 protected:
  static TransProb _instance;
};

#endif // _TRANSPROB_H
