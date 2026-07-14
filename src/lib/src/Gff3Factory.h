//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GFF3FACTORY_H
#define GFF3FACTORY_H
#include "DataSetAbstractFactory.h"
class DataSetAbstractProduct;
#if __GNUC__ < 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ < 3 )
 class ifstream;
#endif

/**
 * <b>Responsibilities:</b><br>
 * instanciation of DataSetAbstractFactory, parses sequence similarity data<br>
 * and stores in the approriate datastructure<br>
 * <br>
 * <b>Collaborators:</b><br>
 * NapDataSet<br>
 * @stereotype factory 
 */
class Gff3Factory : public DataSetAbstractFactory {
public:

  /**
   * format for homology data:
   * 5' 3' relative 5' relative 3' "accession id" "% identiy" "% similarity"
   * 24837 25289 48 500 AU088275 83.12 91.98
   * 
   */
  virtual DataSetAbstractProduct * createDataSetAbstractProduct(const string&, const string&, const string&, const std::string&);

private:
  static DataSetAbstractProduct* gffHomologyRead(std::ifstream&,const string&, const std::string&);

    /**
     * @link
     * @shapeType PatternLink
     * @pattern Abstract Factory
     * @supplierRole Product 
     */
    /*# NapDataSet _genePredDataSet; */
};

#endif //GFF3FACTORY_H
