//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEPREDDATASETCONCRETEFACTORY_H
#define GENEPREDDATASETCONCRETEFACTORY_H
#include "DataSetAbstractFactory.h"
class DataSetAbstractProduct;
class DnaStr;

/**
 * <b>Responsibilities:</b><br>
 * Instance of DataSetAbstractFactory<br>
 * parses a flat file containing one gene prediction data set for <br>
 * a given genomic sequence<br>
 * <br>
 * <b>Collaborators:</b><br>
 * GenePredDataSet, GenePrediction, GeneModelType
 * @stereotype factory 
 */
class GenePredDataSetConcreteFactory : public DataSetAbstractFactory {
public:

  /**
   * data fromat:
   * "gene model number" "exon model" end5 end3
   * example
   * 1 single 7623 7958
   *
   * exon model types are initial, single, terminal, internal
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
    /*# GenePredDataSet _genePredDataSet; */
};
#endif //GENEPREDSETCONCRETEFACTORY_H
