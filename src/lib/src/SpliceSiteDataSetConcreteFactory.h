//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef SPLICESITEDATASETCONCRETEFACTORY_H
#define SPLICESITEDATASETCONCRETEFACTORY_H
#include "DataSetAbstractFactory.h"
class DataSetAbstractProduct;

/**
 * <b>Responsiblities:</b><br>
 * Instanciates the <i>DataSetAbstractFactory</i>, parses and reads in splice site <br>
 * predictions from a flat file and stores in <i>SpliceSiteDataSet</i> object<br>
 * <br>
 * <b>Collaborators:</b><br>
 * SpliceSiteDataSet, SpliceSite
 * @stereotype factory 
 */
class SpliceSiteDataSetConcreteFactory : public DataSetAbstractFactory {
public:

  /**
   * data format:
   * 5' 3' type 
   * 5'>3' --> neg. strnd
   * example:
   * 3200 3219 acceptor
   */
  virtual DataSetAbstractProduct * createDataSetAbstractProduct(const string&, const string&, const string&, const std::string&);

private:

    /**
     * @link
     * @shapeType PatternLink
     * @pattern Abstract Factory
     * @supplierRole Product 
     */
    /*# SpliceSiteDataSet _spliceSiteDataSet; */
};

#endif //SPLICESITEDATASETCONCRETEFACTORY_H
