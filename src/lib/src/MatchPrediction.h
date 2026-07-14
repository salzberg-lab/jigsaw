//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef MATCHPREDICTION_H
#define MATCHPREDICTION_H

#include "GenePrediction.h"
#include <iostream>
class GeneModelType;


/**
 * <b>Responsibilities:</b><br>
 * A further specialization of a gene prediction, based specifically on<br>
 * a sequence similarity match. Contains sequence similarity information,<br>
 * including % identity, % similarity and the local alignment coordinates in <br>
 * the sequence.<br> 
 */
class MatchPrediction : public GenePrediction {
public:
  /**
   * @param const SeqLoc& --> sequence location to genomic sequence
   * @param AnnotatinScore sc --> a confidence score
   * @param int --> associate a user defined model id or -1
   * @param const string& --> string identify the source (e.g. an  accession id)
   * @param double --> % identity of the match
   * @param double --> % similarity of the match
   * @param const SeqLoc& --> the local sequence location on the matching genomic sequence
   * (not accessible at the moment)
   * @param const GeneModelType& --> an exon type (e.g. initial, internal, single, terminal
   */
  MatchPrediction(const SeqLoc& sl, AnnotationScore sc,  int modelId, const string& accession, 
		  double pcntId, double pcntSim, double lscore, const SeqLoc& msl, const GeneModelType& gm) : 
    GenePrediction(sl,sc,modelId,accession,gm), 
		_pcntId(pcntId), _pcntSim(pcntSim), _lscore(lscore), _sl(msl) { }

  virtual ~MatchPrediction() { }

  /** 
   * get the percent identity of the match
   */
  double getPcntId() const { return _pcntId; }

  /**
   * get the percent similarity of the match
   */
  double getPcntSim() const { return _pcntSim; }

  /**
   * get the alignment score
   */
  double getAlignScore() const { return _lscore; }

  /** 
   * standard print
   */
  void print(std::ostream&,const std::string&) const;

private:
  double _pcntId, _pcntSim, _lscore;
  SeqLoc _sl;

};

#endif // SEQMATCH_H
