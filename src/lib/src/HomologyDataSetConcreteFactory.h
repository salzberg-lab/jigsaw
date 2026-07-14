//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef HOMOLOGYDATASETCONCRETEFACTORY_H
#define HOMOLOGYDATASETCONCRETEFACTORY_H
#include "DataSetAbstractFactory.h"
class DataSetAbstractProduct;
class DnaStr;

/**
 * <b>Responsibilities:</b><br>
 * instanciation of DataSetAbstractFactory, parses sequence similarity data<br>
 * and stores in the approriate datastructure<br>
 * <br>
 * <b>Collaborators:</b><br>
 * HomologyDataSet<br>
 * @stereotype factory 
 */
class HomologyDataSetConcreteFactory : public DataSetAbstractFactory {
public:

  /**
   * format for homology data:
   * 5' 3' relative 5' relative 3' "accession id" "% identiy" "% similarity"
   * 24837 25289 48 500 AU088275 83.12 91.98
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
    /*# HomologyDataSet _genePredDataSet; */
};

#endif //HOMOLOGYSETCONCRETEFACTORY_H
