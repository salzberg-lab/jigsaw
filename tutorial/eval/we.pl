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
my %truegenes = loadTrueGenes2($g);

my $fileid = shift( @ARGV );
my $filt = shift( @ARGV );
doEval( $root,$fileid, \%truegenes, $filt ) ;
exit 0;

################

sub doEval {
  my ($root,$fileid,$truegenes,$filt) = @_;
  #my @filelst = @$flstRef;
  my %tg = %$truegenes;
  my $ntg = 0;
  my $we = 0;
  foreach my $seqname (keys %tg) {
    my $tgHelp = $tg{$seqname};
    my %justOnce;
    my %intronOnce;
    my %trueExon;
    foreach my $geneid (keys %$tgHelp) {
      my $tgene = $tgHelp->{$geneid};
      my @vals=split(/ /,$tgene);
      for(my $cnt=1; $cnt <= $#vals; $cnt++) {
	my $exon = "$vals[$cnt] $vals[$cnt+1]";
	$trueExon{$exon}=1;
      }
    }
    foreach my $geneid (keys %$tgHelp) {
      my $tgene = $tgHelp->{$geneid};
      if( $justOnce{$tgene} ) {
	next;
      }
      $justOnce{$tgene} = 1;
      my %exonTrack;
      my %geneTrack;
      #my @tgenes=split(/:/,$truegenes->{$seqname});
      my @tgenes;
      push(@tgenes,$tgene);
      my $file = "$root/$seqname/$seqname.$fileid";
      if ( ! (-e $file) ) {
	$file = "$root/$seqname/${seqname}_$fileid";
      }
      my $predgenes=loadGff($file,$filt);
      my @pgenes=split(/:/,$predgenes);
      $we += wrongExons($seqname,$fileid,\@pgenes,$tgene,\%geneTrack,\%intronOnce,\%trueExon);
      $ntg += $#tgenes +1;
    }
  }
  #print "WE: $we\n";
}

sub wrongExons {
  my ($asmbl, $ptype, $pgenes, $tg, $geneTrack, $intronOnce, $trueExon) = @_;
  my $we =0;
  my @mytvals = split(/:/,$tg);
  for(my $myiter = 0; $myiter < $#mytvals + 1; $myiter++) {
    my @tintrons=split(/ /,$mytvals[$myiter]);
    my $tstrnd = $tintrons[0];
    for(my $introni = 2; $introni < $#tintrons; $introni+=2) {
      my ($i5,$i3) = ($tintrons[$introni],$tintrons[$introni+1]);
      my $intron = "$i5 $i3";
      if( $intronOnce->{$intron} ) {
	next;
      }
      $intronOnce->{$intron}=1;
      for(my $i=0;$i<@$pgenes;$i++) {
	my @exons = split(/ /,$pgenes->[$i]);
	my $pstrnd = $exons[0];
	if( $pstrnd eq $tstrnd ) {
	  for(my $ei = 1; $ei <= $#exons-1; $ei+=2) {
	    my ($end5,$end3) = ($exons[$ei],$exons[$ei+1]);
	    #print "$i5 $i3 $end5 $end3\n";
	    if( $end5 > $i5 && $end5 < $i3 && $end3 > $i5 && $end3 < $i3 && !$trueExon->{"$end5 $end3"} ) {
	      print "we: $end5 $end3\n";
	      $we++;
	    }
	  }
	}
      }
    }
  }
  return $we;
}
