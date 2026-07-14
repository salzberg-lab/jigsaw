#!/usr/bin/perl

use strict;

my $score = 1.0;
my @genes;

while(my $line = <STDIN>) {
  my @vals = split(/\s+/,$line);
  my $id = $vals[0];
  my ($idx1,$idx2)=(2,3);
  if( $#vals == 2 ) {
     ($idx1,$idx2)=(1,2);
  } else {
	  $id .= " $vals[1]";
  }
  my $end5 = $vals[$idx1];
  my $end3 = $vals[$idx2];
  my $strnd = "+";
  if( $end5 > $end3 ) {
    ($end5,$end3)=($end3,$end5);
    $strnd = "-";
  }
  my $type="exon";
  my $score=1.0;
  my $frame=1;
  push (@genes, "seqid\tgp\t$type\t$end5\t$end3\t$score\t$strnd\t$frame\t$id\n");
}

foreach my $line (sort { (split(/\t/,$a))[2] <=> (split(/\t/,$b))[2] } @genes) {
  print $line;
}
