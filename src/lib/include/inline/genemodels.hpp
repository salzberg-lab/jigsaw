//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
// end5 should be end3 here darn it
bool genemodels::isContiguous(int end5, dsu::Strand_t s5, 
			      int end3, dsu::Strand_t s3) {
  return (end5+1)==end3 && s5==s3;
}
