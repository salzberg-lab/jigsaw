#!/usr/bin/perl

# program takes all of the computational gene prediction methods
# including combiners, excluding sequence homology
# tracks which methods got which gene predictions right.
# operates only on exons and whole genes

use strict;
use CombUtils;


################
## main        #
################

# f=rootdirectory for $gp file
# g=true gene file
# filetype=geneprediction file name
my $root = shift( @ARGV );
my $g = shift( @ARGV );
my $doGenes = shift( @ARGV );
my $fileid = shift( @ARGV );
my $filt=shift(@ARGV);
my $genesorexons = "exons";
if( $doGenes ) {
  $genesorexons = "genes";
}

#print "Compute number of completely missed $genesorexons\n";
    my %truegenes = loadTrueGenes2($g);

    doEval( $root, $fileid, \%truegenes, $doGenes, $filt ) ;
    exit 0;

################

sub doEval {
  my ($root,$fileid,$truegenes,$doGene,$filt) = @_;
  #my @filelst = @$flstRef;
  my %tg = %$truegenes;
  my $ntg = 0;
  foreach my $seqname (keys %tg) {
    my $tgHelp = $tg{$seqname};
    foreach my $geneid (keys %$tgHelp) {
      my $tgene = $tgHelp->{$geneid};
      #print "hmm [$tgene]\n";
      #exit(1);
      my %exonTrack;
      my %geneTrack;
      #my @tgenes=split(/:/,$truegenes->{$seqname});
      my @tgenes;
      push (@tgenes,$tgene);
      my @texf = compTextf(@tgenes);
      my $file = "$root/$seqname/$seqname.$fileid";
      #print "check: $file\n";
      if ( ! (-e $file) ) {
	$file="$root/$seqname/${seqname}_$fileid";
      }
      if ( ! (-e $file)) {
	print "could not open $file unable to continue\n";
	exit(1);
      }
      my $predgenes=loadGff($file,$filt);
      my @pgenes=split(/:/,$predgenes);
      if ( $doGene ) {
	compCgenes($seqname,$fileid,\@pgenes,$tgene,\%geneTrack);
      } else {
	compCexons($seqname,$fileid,\@tgenes,\@pgenes,\%exonTrack);
      }
      if ( $doGene ) {
	foreach my $key (keys %geneTrack) {
	  print "$key:$geneTrack{$key}\n";
	}
      } else {
	foreach my $key (keys %exonTrack) {
	  print "$key:$exonTrack{$key}\n";
	}
      }
      $ntg += $#tgenes +1;
    }
  }
  #print "$ntg\n";
}

sub compTextf {
  my (@tgenes) = @_;
    my @texf;
    $texf[0]="";$texf[1]="";$texf[2]="";$texf[3]="";$texf[4]="";$texf[5]="";
    my $atex="";
    my $apex="";
    my $cef=0;
    for(my $i=0;$i<=$#tgenes;$i++) {
	my $frame=0;
	my @e=split(/\s+/,$tgenes[$i]);
	for(my $j=1;$j<=$#e;$j+=2) {
	    if($e[0] eq "+") {
		$texf[$frame].=" ".$e[$j]." ".$e[$j+1]." ";
	    }
	    else {
		$texf[3+$frame].=" ".$e[$j]." ".$e[$j+1]." ";
	    }
	    $atex.=" ".$e[$j]." ".$e[$j+1]." ";
	    $frame=($frame + $e[$j+1] - $e[$j] + 1) % 3;
	}
    }
  return @texf;
}

sub cntPrcds {
  my $begseq = shift @_;
  my $endseq = shift @_;
  my (@b) = @_;

    # $prcds = no. of predicted coding base pairs
    my $prcds=0;
    my $i=0;
    while($i<=$#b) {
	if($b[$i] =~ /^\d+$/) { # $b[$i] is a number
	}
	else {
	    $i++;
	}
	my $l1=$b[$i];$i++;
	if($l1 < $begseq) { $l1=$begseq;}
	if($l1 > $endseq) { $l1=$endseq+1;}
	my $l2=$b[$i];$i++;
	if($l2 < $begseq) { $l2=$begseq-1;}
	if($l2 > $endseq) { $l2=$endseq;}
	$prcds+=$l2-$l1+1;
	#print "$prcds\n";
	#if($prcds < 0) {
	    #exit (0);
	#}
    }
  return $prcds;
}

sub compCexons {
  my ($gid,$pid,$texfRef,$pgenesRef,$exonTrackRef) = @_;
  my @tgenes = @$texfRef;
  my @pgenes = @$pgenesRef;
  #my %exonTrack = %$exonTrackRef;
  for(my $iter=0;$iter<=$#tgenes;$iter++) {
    my @texons = split(/\s+/,$tgenes[$iter]);
    my $tstrnd = $texons[0];
    for(my $eiter = 1; $eiter < $#texons; $eiter+=2) {
      my $tend5 = $texons[$eiter];
      my $tend3 = $texons[$eiter+1];
      my $isLap = 0;
      for(my $i=0;$i<=$#pgenes;$i++) {
	#print "wait: [$pgenes[$i]]\n";
	my @pexons = split(/\s+/,$pgenes[$i]);
	my $pstrnd = $pexons[0];
	if( $pstrnd eq $tstrnd ) {
	  for (my $pe = 1; $pe < $#pexons; $pe+=2) {
	    my $pend5 = $pexons[$pe];
	    my $pend3 = $pexons[$pe+1];
	    my $tl = "$tstrnd $tend5 $tend3";
	    #print "huh? $pe [$pend5 $pend3] $tl\n";
	    if (isOverlap($pend5,$pend3,$tl) ) {
	      $isLap = 1;
	      last;
	    }
	  }
	}
	if( $isLap ) {
	  last;
	}
      }
      if( !$isLap ) {
	my $id = "$gid,$tend5,$tend3";
	$exonTrackRef->{$id} .= $pid . " ";
      }
    }
  }
}

sub isIn {
  my ($val,$pair) = @_;
  my @vals = split(/ /,$pair);
  if( $val >= $vals[1] && $val <= $vals[$#vals] ) {
    #print "true? $val, $vals[1], $vals[$#vals]\n";
    return 1;
  }
  return 0;
}

sub isOverlap {
  my ($left,$right,$pair) = @_;
  my @vals = split(/ /,$pair);
  if( ($left < $vals[1] && $right > $vals[$#vals]) ||
      isIn($left,$pair) || isIn($right,$pair) ) {
    return 1;
  }
  return 0;
}


sub compCgenes {
  my ($asmbl, $ptype, $pgenes, $tg, $geneTrack) = @_;
  my @mytvals = split(/:/,$tg);
  #print "FHDDSFDS: $ptype\n";
  #print "@mytvals\n";
  for(my $myiter = 0; $myiter < $#mytvals + 1; $myiter++) {
    my $found = 0;
    #print "true: vals: [$mytvals[$myiter]]\n";
    my $tstrnd = "-";
    if( $mytvals[$myiter] =~ /\+/ ) {
      $tstrnd = "+";
    }
    for(my $i=0;$i<@$pgenes;$i++) {
      #print "$pgenes->[$i]\n";
      if( ($pgenes->[$i] !~ /\+/ && $tstrnd eq "+" )  ||
	  ($pgenes->[$i] !~ /\-/ && $tstrnd eq "-") ){
	next;
      }
      my @val = split(/ /,$pgenes->[$i]);
      my $pstrnd = $val[0];
      my $pend5 = $val[1];
      my $pend3 = $val[$#val];
      #print "pred: vals: [@val]\n";
      #print "[$tstrnd] $pstrnd\n";
      #print "WHATUP: $pend5, $pend3\n";
      if($tstrnd eq $pstrnd && $mytvals[$myiter]eq$pgenes->[$i]) {
	#print "pg: [$pgenes[$i]]\ntg: [$truegenes{$seqname}]\n";
	#print "pg: [$pgenes->[$i]]\ntg: [$mytvals[$myiter]]\n";
	#print "overlap!!!!!\n";
	#my $id = "$asmbl $mytvals[$myiter]";
	#$geneTrack->{$id} .= $ptype . " ";
	#$ncorr++;
	$found = 1;
	last;
      } elsif($tstrnd eq $pstrnd && isOverlap($pend5,$pend3,$mytvals[$myiter]) ) {
	#print "pg: [$pgenes->[$i]]\ntg: [$mytvals[$myiter]]\n";
	#print "overlap!!!!!\n";
	#my $id = "$asmbl,$mytvals[$myiter],$pgenes->[$i]";
	#$geneTrack->{$id} .= $ptype . " ";
	#last;
	$found = 1;
	last;
      }
    }
    if( !$found ) {
      my $id = "$asmbl,$mytvals[$myiter]";
      #print "Not found: $ptype\n";
      $geneTrack->{$id} .= $ptype . " ";
      #exit(1);
    }
  }
}
