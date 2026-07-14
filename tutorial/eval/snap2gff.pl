#!/usr/bin/perl

use strict;

my @genes;
while(my $line=<STDIN>) {
  chomp($line);
  if( $line =~ />/ ) {
    next;
  }
  my @vals =split(/\s+/,$line);
  my $end5=$vals[1];
  my $end3=$vals[2];
  my $strnd=$vals[3];
  my $score =$vals[4];
  my $id =$vals[8];
  my @tv=split(/\./,$id);
  my $gnum=$tv[1];
  #if($strnd eq"-") {
    #($end5,$end3)=($end3,$end5);
  #}
  my $exon=$vals[0];
  #print "$gnum $exon $end5 $end3 $score\n";
  my $frame=1;
  push @genes, "seqid\tsnap\t$exon\t$end5\t$end3\t$score\t$strnd\t$frame\t$id\n";
}

foreach my $line (sort { (split(/\t/,$a))[3] <=> (split(/\t/,$b))[3] } @genes) {
	print $line;
}
