#!/usr/bin/perl -w

use strict;
use CombUtils;
use MyUtils;

my $root=$ARGV[0] || die "need root directory\n";
my $trueFile=$ARGV[1] || die "need answer genes\n";
my $predFileId=$ARGV[2] || die "need pred genes\n";
my $filt=$ARGV[3];
if( !defined($filt) || !$filt) {
	$filt="";
}

my %trueGenes = loadTrueGenes2($trueFile);
my %double=();
my @rs = (0,0,0,0,0,0,0,0,0,0,0,0);
my $psum=0;
my $alts=0;
my $loci=0;
foreach my $did (keys %trueGenes) {
  my $fastaFile = "$root/$did/$did.contig";
  if( ! -e $fastaFile && -e "$root/$did/$did.seq" ) {
    $fastaFile = "$root/$did/$did.seq";
  } elsif( ! (-e $fastaFile) && -e "$root/$did/$did.fa")  {
    $fastaFile = "$root/$did/$did.fa";
  }
  my $pip = "tail +2 $fastaFile | chop.pl | wc |";
  open(PIPE,$pip) || die "failed [$pip]\n";
  my $val = <PIPE>;
  chomp($val);
  my $slen=0;
  if ( $val =~ /\s+\d+\s+\d+\s+(\d+)/ ) {
    $slen = $1;
  }
  close(PIPE);
  my $tgHelp = $trueGenes{$did};
  my @tmpgene;
  my %onlyOnce;
  foreach my $geneid (keys %$tgHelp) {
    my $tgenes = $tgHelp->{$geneid};
    if( !$onlyOnce{$tgenes} ) {
      push(@tmpgene,$tgenes);
      $onlyOnce{$tgenes}=1;
    }
  }
  my @tgs=();
  my @sgene = sort { (split(/ /,$a))[1] <=> (split(/ /,$b))[1] } @tmpgene;
  my $iter=0;
  for($iter = 0; $iter <= $#sgene; ) {
    my @curr=split(/ /,$sgene[$iter]);
    my $cset = $sgene[$iter];
    my $nval=$iter+1;
    my $lValid = $nval;
    my $lv=0;
    while($nval <= $#tmpgene) {
      my @nurr=split(/ /,$sgene[$nval]);
      #print "check $sgene[$nval]\n";
      if( $curr[0] eq $nurr[0] && $nurr[1] < $curr[$#curr] ) {
	$cset .= ":$sgene[$nval]";
	my $v=$sgene[$nval];
	#if( $cset =~ /$v/  ) {
	  #print "add: $cset\n";
	  #exit(1);
	#}
	$lv=1;
	$lValid = $nval;
	$nval++;
      } else {
	last;
      }
    }
    if($lv) {
      $lValid++;
    }
    push(@tgs,$cset);
    my $range = ($lValid-1);
    if( $range > $iter ) {
      $alts++;
    }
    #print "genes in same loci: $iter $range\n";
    $iter = $lValid;
  }
  my $predFile="$root/$did/$did.$predFileId";
  if ( !(-e $predFile) ) {
    $predFile="$root/$did/${did}_$predFileId";
  }
  my $pgenes = loadGff($predFile,$filt);
  #print "crapola: $pgenes\n";
  #exit(1);
  print "proc: $predFile\n";
  my %onlyOnce1;
    my %pgeneHash=();
    my %texons=();
    my %tintrons=();
    my %tnuc=();
    my %pexons=();
    my %pintrons=();
    my %noDups=();
  for(my $iter = 0; $iter<=$#tgs; $iter++) {
    my $rghtSide=-1;
    my $leftSide=-1;
    if( $iter > 0 ) {
      my @ttgenes = split(/:/,$tgs[$iter-1]);
      my @tmp = split(/ /,$ttgenes[$#ttgenes]);
      $leftSide = $tmp[$#tmp];
    }
    if( $iter < $#tgs ) {
      my @ttgenes = split(/:/,$tgs[$iter+1]);
      my @tmp = split(/ /,$ttgenes[$#ttgenes]);
      $rghtSide = $tmp[$#tmp];
    }
    my $tga=$tgs[$iter];
    $loci++;
    my @tgenes = split(/:/,$tga);
    my $br = (split(/ /,$tgenes[0]))[1]  - 200;
    my @tmp = split(/ /,$tgenes[$#tgenes]);
    my $er = $tmp[$#tmp]+200;
    if($er > $slen) {
      $er = $slen-1;
    }
    if( $leftSide != -1 && $br < $leftSide ) {
      print "fixing er problem $br $leftSide\n";
      $br=$leftSide+1;
    }
    if( $rghtSide != -1 && $er > $rghtSide ) {
      print "fixing er problem $er $rghtSide\n";
      $er=$rghtSide-1;
    }
    foreach my $tgene (@tgenes) {
      if( $onlyOnce1{$tgene}) {
	print "what the f: $tga\n";
	exit(1);
      }
      if( !$onlyOnce1{$tgene} ) {
	$onlyOnce1{$tgene}=1;
      }
      my ($tg,$tex,$ti,$tn) = cntTrue($tgene,\%texons,\%tintrons,\%tnuc);
      my @results =
	cntOverPred($tgene,$pgenes,$br,$er,\%pexons,\%pintrons,\%noDups,\%tnuc,\%pgeneHash);
      #print "$tg $tgene\n";
      ($results[8],$results[9],$results[10],$results[11]) = ($tg,$tex,$ti,$tn);

      for(my $iter= 0; $iter <= $#results; $iter++) {
	$rs[$iter] += $results[$iter];
	#print "crap:$iter $rs[$iter] $results[$iter]\n";
      }
    }
  }
}

my ($numPredGenes,$numPredExons,$numPredIntrons,$numPredNuc) = ($rs[0],$rs[1],$rs[2],$rs[3]);
my ($numCorrGenes,$numCorrExons,$numCorrIntrons,$numCorrNuc) = ($rs[4],$rs[5],$rs[6],$rs[7]);
my ($totGenes,$totExons,$totIntrons,$totNuc)                 = ($rs[8],$rs[9],$rs[10],$rs[11]);
my $sng = $numCorrGenes/$totGenes;
my $spg = $numCorrGenes/$numPredGenes;

my $sne = $numCorrExons/$totExons;
my $spe = $numCorrExons/$numPredExons;

my $sni = $numCorrIntrons/$totIntrons unless ($totIntrons < 1);
my $spi = $numCorrIntrons/$numPredIntrons unless ($totIntrons < 1);

my $snn = $numCorrNuc/$totNuc;
my $spn = $numCorrNuc/$numPredNuc;
print "total predicted number of genes: $psum, alts: $alts $loci\n";
print "raw $numPredGenes $numPredExons $numPredIntrons $numPredNuc ";
print "$numCorrGenes $numCorrExons $numCorrIntrons $numCorrNuc ";
print "$totGenes $totExons $totIntrons $totNuc\n";
print "sng: $sng spg: $spg ";
print "sne: $sne spe: $spe ";
print "sni: $sni spi: $spi ";
print "snn: $snn spn: $spn\n";

sub cntTrue {
  my ($tgene,$texRef,$tinRef,$tnucRef)=@_;
  my @tcoords = split(/ /,$tgene);
  my ($geneCnt,$exCnt,$inCnt,$nucCnt) = (1,0,0,0);
  for(my $iter =1; $iter <= $#tcoords; $iter+=2) {
    my $exon = "$tcoords[$iter] $tcoords[$iter+1]";
    if( !$texRef->{$exon} ) {
      $exCnt++;
      $texRef->{$exon} = 1;
    }
    for(my $j = $tcoords[$iter]; $j <= $tcoords[$iter+1]; $j++) {
      if( !$tnucRef->{$j} ) {
	$nucCnt++;
	$tnucRef->{$j}=1;
      }
    }
    if( $iter + 2 < $#tcoords ) {
      my $intr = "$tcoords[$iter+1] $tcoords[$iter+2]";
      if( !$tinRef->{$intr} ) {
	$inCnt++;
	$tinRef->{$intr} = 1;
      }
    }
  }
  #print "huh? $geneCnt $exCnt $inCnt $nucCnt\n";
  return ($geneCnt,$exCnt,$inCnt,$nucCnt);
}

sub cntOverPred {
  my ($tgene,$pgeneStr,$br,$er,$pexonsRef,$pintronsRef,$noDupsRef,$trueNuc,$pgHash) = @_;
  my @pgenes = split(/:/,$pgeneStr);
  my ($numPredGenes,$numPredExons,$numPredIntrons,$numPredNuc) = (0,0,0,0);
  my ($numCorrGenes,$numCorrExons,$numCorrIntrons,$numCorrNuc) = (0,0,0,0);
  my ($totGenes,$totExons,$totIntrons,$totNuc)=(0,0,0,0);
  my $prev=-1;
  my @tcoords = split(/ /,$tgene);
  my $tstrnd = shift(@tcoords);
  my @tsort = sort {$a <=> $b} @tcoords;
  #print "tg: $tgene\n";
  foreach my $pgene ( @pgenes ) {
    if ( !$pgene ) {
      next;
    }
    my @pcoords = split(/ /,$pgene);
    my $pstrnd = shift(@pcoords);
    my @psort = sort {$a <=> $b} @pcoords;
    my ($pb,$pe) = ($psort[0],$psort[$#psort]);
    #print "pg: $pgene\n";
    #print "huh? [$pb] [$pe] [$tb] [$te] [$er] [$br]\n";
    if ((($pb >= $br && $pb <= $er) ||
	 ($pe >= $br && $pe <= $er) ||
	 ($pb <= $br && $pe >= $er)) ) {
      if ( !$pgHash->{$pgene} ) {
	$numPredGenes++;
	#print "crap! [$pgene]\n";
	$pgHash->{$pgene} = 1;
      }
      #print "shot: $pgene\n";
      if ( $pgene eq $tgene ) {
	#print "humm:\n$pgene\n$tgene\n";
	$numCorrGenes++;
      }
      my ($ce,$ci,$cn,$peCnt,$piCnt,$pn) = countCorrExons($pstrnd,$tstrnd,\@tsort,\@psort, $pexonsRef,
							  $pintronsRef,$noDupsRef,$trueNuc);
      $numCorrExons += $ce;
      $numCorrIntrons += $ci;
      $numCorrNuc += $cn;
      #push(@overlaps,\@psort);

      $numPredExons += $peCnt;
      $numPredIntrons += $piCnt;
      $numPredNuc += $pn;
    }
  }
  #print "woa: $numPredGenes $numPredExons $numPredIntrons $numPredNuc\n";
  return ($numPredGenes,$numPredExons,$numPredIntrons,$numPredNuc,
	  $numCorrGenes,$numCorrExons,$numCorrIntrons,$numCorrNuc,
	  $totGenes,$totExons,$totIntrons,$totNuc);
}

sub remOverlapsPred {
  my ($ovRef) = @_;
  my @overlaps = @$ovRef;
  my $remOv=0;
  my ($min,$max)=(-1,-1);
  for(my $iter = 0; $iter <= $#overlaps; $iter++) {
    my $tref = $overlaps[$iter];
    my @p1 = @$tref;
    my ($l1,$r1) = ($p1[0],$p1[$#p1]);
    #print "hm: [$l1 $r1] [$min $max]\n";
    if( $min != -1 && $max != -1) {
      if( ($l1 >= $min && $l1 <= $max) ||
	  ($r1 >= $min && $r1 <= $max) ||
	  ($min >= $l1 && $min <= $r1) ||
	  ($max >= $l1 && $max <= $r1))
	  {
	    $remOv++;
	    if( ($max-$min) < ($r1-$l1) ) {
	      ($min,$max)=($l1,$r1);
	    }
	  }
    } else {
      ($min,$max) = ($p1[0],$p1[$#p1]);
    }
  }
  return $remOv;
}

sub countCorrExons {
  my ($pstrnd,$tstrnd,$tex,$pex,$pexons,$pintrons,$noDupNuc,$truePrc) = @_;
  my @texons = @$tex;
  my @pexons = @$pex;
  my ($predNucCnt,$corrNucCnt)=(0,0);
  my $pexCnt = 0;
  my $pinCnt = 0;
  my ($cexCnt,$cinCnt) = (0,0);
  for(my $piter = 0; $piter < $#pexons; $piter+=2 ) {
    my $pexon = "$pexons[$piter] $pexons[$piter+1]";
    #print "h: [$pexon]\n";
    my $pintron = "";
    if( $piter+2 < $#pexons ) {
      $pintron = "$pexons[$piter+1] $pexons[$piter+2]";
    }
    for(my $nuc = $pexons[$piter]; $nuc <= $pexons[$piter+1]; $nuc++) {
      if( !$noDupNuc->{$nuc} ) {
	if( $truePrc->{$nuc} && $pstrnd eq $tstrnd) {
	  $corrNucCnt++;
	  $noDupNuc->{$nuc}=2;
	} else {
	  $noDupNuc->{$nuc}=1;
	}
	$predNucCnt++;
      } elsif ($noDupNuc->{$nuc} == 1 && $truePrc->{$nuc} && $tstrnd eq $pstrnd ) {
	$corrNucCnt++;
	$noDupNuc->{$nuc}=2;
      }
    }
    for(my $titer = 0; $titer < $#texons; $titer+=2 ) {
      my $texon = "$texons[$titer] $texons[$titer+1]";
      my $tintron = "";
      if( $titer+2 < $#texons ) {
	$tintron = "$texons[$titer+1] $texons[$titer+2]";
      }
      #print "h: p:$pintron t:$tintron\n";
      if( !$pexons->{$pexon} ) {
	if( $texon eq $pexon && $tstrnd eq $pstrnd ) {
	  $cexCnt++;
	  $pexons->{$pexon} = 2;
	} else {
	  $pexons->{$pexon} = 1;
	}
	$pexCnt++;
      } elsif( $pexons->{$pexon} == 1 && $texon eq $pexon && $tstrnd eq $pstrnd ) {
	$cexCnt++;
	$pexons->{$pexon} = 2;
      }
      if( !$pintrons->{$pintron} ) {
	if( $tintron && $pintron && $tintron eq $pintron && $tstrnd eq $pstrnd ) {
	  $cinCnt++;
	  $pintrons->{$pintron} = 2;
	} else {
	  $pintrons->{$pintron} = 1;
	}
	$pinCnt++;
      } elsif( $pintrons->{$pintron}==1 &&  $tintron && $pintron && $tintron eq $pintron && $tstrnd eq $pstrnd ) {
	$cinCnt++;
	$pintrons->{$pintron} = 2;
      }
    }
  }
  #print "weird: $cexCnt $cinCnt, $pexCnt $pinCnt\n";
  #print "@pexons\n@texons\n";
  return ($cexCnt,$cinCnt,$corrNucCnt,$pexCnt,$pinCnt,$predNucCnt);
}

sub countCorrNuc {
  my ($tex,$apex,$br,$er,$pnucs) = @_;
  my $nsub=0;
  my @texons = @$tex;
  my @paexons = @$apex;
  my @pexons=();
  foreach my $pref (@paexons) {
    @pexons = (@pexons,@$pref);
  }
  my $corrNuc=0;
  my %crds;
  foreach my $crd (@texons) {
    $crds{$crd}=1;
  }
  foreach my $crd (@pexons) {
    if( $crd < $br ) {
      $crd = $br;
    }
    if($crd > $er ) {
      $crd = $er;
    }
    $crds{$crd}=1;
  }
  my @scrds = sort { $a <=> $b } keys %crds;
  for(my $siter = 0; $siter < $#scrds; $siter++) {
    my ($from,$to) = ($scrds[$siter],$scrds[$siter+1]);
    #print "check: $from $to\n";
    my ($match,$isUtr) = isSameOverlap($from,$to,$tex,$apex,\@pexons,$br,$er);
    if( $match ) {
      my $tv = ($to-$from)+1;
      #print "match: $tv\n";
      for(my $cnt = $from; $cnt <= $to; $cnt++) {
	if( !$pnucs->{$cnt} ) {
	  $corrNuc++;
	  $pnucs->{$cnt} = 1;
	} else {
	  $nsub++;
	}
      }
      #$corrNuc+= $tv;
    }
  }
  return ($corrNuc,$nsub);
}

sub isSameOverlap {
  my ($from,$to,$tex,$pex,$pref,$br,$er) = @_;
  my @texons = @$tex;
  my $prev=-1;
  my %plabel;
  my @pexons = @$pref;
  #print "what the hell: [@pexons]\n";
  my ($ptype,$ttype)=("","");
  foreach my $pref (@$pex) {
    my @pex = @$pref;
    if( $prev != -1 ) {
      $plabel{"$prev $pex[0]"} = "IGR";
    }
    for(my $iter = 0; $iter < $#pex; $iter++) {
      my ($e5,$e3)=($pex[$iter],$pex[$iter+1]);
      if( $e5 < $br ) {
	$e5 = $br;
      }
      if($e3 > $er ) {
	$e3 = $er;
      }
      if( ($iter % 2) == 0 || $iter==0 ) {
	$plabel{"$e5 $e3"} = "exon";
      } else {
	$plabel{"$e5 $e3"} = "intron";
      }
      my $tid="$e5 $e3";
      #print "what the f: $tid $plabel{$tid}\n";
    }
    $prev=$pex[$#pex];
  }
  for(my $iter = 0; $iter < $#texons; $iter++) {
    my ($e5,$e3)=($texons[$iter],$texons[$iter+1]);
    if( $from >= $e5 && $from <= $e3 &&
	$to >= $e5 && $to <= $e3 ) {
      if( ($iter % 2) == 0 || $iter==0 ) {
	$ttype = "exon";
      } else {
	$ttype = "intron";
      }
      #print "$e5 $e3 $ttype\n";
      last;
    }
  }
  for(my $iter = 0; $iter < $#pexons; $iter++) {
    my ($e5,$e3)=($pexons[$iter],$pexons[$iter+1]);
    #print "huh? $e5 $e3 $from $to\n";
    if( $from >= $e5 && $from <= $e3 &&
	$to >= $e5 && $to <= $e3 ) {
      $ptype = $plabel{"$e5 $e3"};
      #print "here: we are [$ptype], $e5 $e3\n";
      last;
    }
  }
  my $val = 0;
  my $isUtr=0;
  #print "[$ptype] [$ttype]\n";
  #exit(1);
  if ( !$ptype  || !$ttype ) {
    #print "must be utr: [$from][$to][$ptype][$ttype]\n";
    #print "@texons\n@pexons\n";
    $isUtr=1;
  } elsif ( $ptype eq $ttype && $ptype eq "exon" ) {
    $val = 1;
  }
  return ($val,$isUtr);
}
