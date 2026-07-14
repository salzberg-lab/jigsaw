TraceBackHelp::tback
TraceBackHelp::getPred(int frame, dsu::Strand_t strnd, dsu::side_t side, int type) const {
  const Options& opt = *Options::getInstance();
  const Compare* cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    cmp = GtCompare::getInstance();
  assert(type == getBest(strnd,frame,side,type).first);
  const AnnotationItem* pred = getBest(strnd,frame,side,type).second.first;
  if( !pred ) {
    return tback(NULL,strnd_frm_t(dsu::eEither,-1),eSide,type);
  }
  dsu::side_t nside = side==dsu::eLeft?dsu::eRght:dsu::eLeft;
  dsu::Strand_t nstrnd = dsu::eEither;
  int nfrm = frame==0?0:3-frame, ntype = -1;
  AnnotationScore bscore=AnnotationScore::badVal();
  for(unsigned ltype = 0; ltype < 6; ++ltype) {
    for(unsigned i = 0; i < 2; ++i) {

  bool someValid=false;
  for(unsigned ltype = 0; ltype < 6; ++ltype) {
    for(unsigned i = 0; i < 2; ++i) {
      dsu::Strand_t st = (dsu::Strand_t)i;
      if(type == Igr::getInstance().getScoreIdx()) {
	lside = eRght;
	if(pred->getBest(st,0,lside,ltype).second.second.isValid()) {
	  if(cmp->compare(pred->getBest(dsu::ePos,0,lside,ltype).second.second,bscore)) {
	    nstrnd = dsu::ePos;
	    bscore = pred->getBest(dsu::ePos,0,lside,ltype).second.second;
	    ntype = ltype;
	  }
	  nfrm = 0;
	  someValid=true;
	}
      }
      if( type == Intron::getInstance().getScoreIdx()) {
	lside = eRght;
	nstrnd = strnd;
	nfrm = frame > 0 ? 3-(frame % 3) : frame;
	someValid=true;
  } else if( !someValid ) {
    lside = eLeft;
    assert(side==eRght);
    AnnotationScore bscore = AnnotationScore::badVal();
    int tfrm=-1;
    for(unsigned ltype = 0; ltype < 6; ++ltype) {
    if(type==Terminal::getInstance().getScoreIdx()) {
      if(strnd==dsu::ePos) {
	tfrm = (getEnd3()-pred->getEnd5()+1) % 3;
      } else {
	tfrm=0;
      }
    } else if(type==Initial::getInstance().getScoreIdx()) {
      if(strnd==dsu::ePos) {
	tfrm = 0;
      } else {
	tfrm = (getEnd3()-pred->getEnd5()+1) % 3;
      }
    } else if(type==Single::getInstance().getScoreIdx()) {
      tfrm=0;
    } else if(type==Internal::getInstance().getScoreIdx()) {
      int len = getEnd3()-pred->getEnd5()+1;
      tfrm = abs(len-frame) % 3;
    }
    if(pred->getBest(strnd,tfrm,lside,ltype).second.second.isValid()) {
      if(cmp->compare(pred->getBest(strnd,tfrm,lside,ltype).second.second,bscore)) {
	nstrnd = strnd;
	someValid=true;
	assert(tfrm!=-1);
	nfrm=tfrm;
      }
    }
    }
  }
  assert(someValid);
  assert(lside != eSide);
  assert(pred->getBest(nstrnd,nfrm,lside,ntype).second.second.isValid());
  //int ntype = pred->getBest(nstrnd,nfrm,lside,ltype).first;
  return tback(pred,strnd_frm_t(nstrnd,nfrm),lside,ntype);
}
