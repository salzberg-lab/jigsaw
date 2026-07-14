//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEMARKHMMFACTORY_H
#define GENEMARKHMMFACTORY_H
#include "DataSetAbstractFactory.h"
class DataSetAbstractProduct;
class ifstream;

/**
 * <b>Responsibilities:</b><br>
 * instanciation of DataSetAbstractFactory, parses sequence similarity data<br>
 * and stores in the approriate datastructure<br>
 * <br>
 * <b>Collaborators:</b><br>
 * NapDataSet<br>
 * @stereotype factory 
 */
class GenemarkHmmFactory : public DataSetAbstractFactory {
public:

  /**
   * reads GenemarkHmm formatted files (geneprediction type)
   * 
   */
  virtual DataSetAbstractProduct * createDataSetAbstractProduct(const string&, const string&, const string&, const std::string&);

private:

    /**
     * @link
     * @shapeType PatternLink
     * @pattern Abstract Factory
     * @supplierRole Product 
     */
    /*# NapDataSet _genePredDataSet; */
};

#endif //GENEMARKHMMFACTORY_H
