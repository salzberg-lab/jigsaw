//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef INTERNAL_H
#define INTERNAL_H

#include "GeneModelType.h"

/**
 * <b>Responsibilities:</b><br>
 * instantiation of <i>GeneModelType</i>, defines behavior<br>
 * of an "internal" exon.<br> 
 */
class Internal : public GeneModelType {
protected:

    Internal();

public:
	void print(std::ostream&,const string&) const;
	const string& getTypeStr() const;
	const int getScoreIdx() const { return 1; }

	static const Internal& getInstance() { return _internal; }
private:
	static Internal _internal;

    /**
     * @link
     * @shapeType PatternLink
     * @pattern Singleton
     * @supplierRole Singleton factory 
     */
    /*# Internal _internal1; */
};

#endif // INTERNAL_H
