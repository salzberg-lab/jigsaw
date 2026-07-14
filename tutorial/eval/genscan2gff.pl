#!/usr/bin/perl

use strict;

my $score = 1.0;
my @genes;

while(my $line = <STDIN>) {
  my @vals = split(/\s+/,$line);
    my ($id,$type,$strnd,$end5,$end3,$frame,$score) = (int($vals[0]),$vals[1],$vals[2],$vals[3],$vals[4],$vals[6],$vals[11]);
    if( $type eq "Init" ) {
      $type = "Initial";
    }
    elsif( $type eq "Intr" ) {
      $type = "Internal";
    }
    elsif( $type eq "Term" ) {
      $type = "Terminal";
    }
    elsif( $type =~ /S/ ) {
      $type = "Single";
    } else {
      next;
    }
    if( $end5 > $end3) {
      ($end5,$end3)=($end3,$end5);
    }
    push (@genes, "seqid\tgenscan+\t$type\t$end5\t$end3\t$score\t$strnd\t$frame\t$id\n");
  }

foreach my $line (sort { (split(/\t/,$a))[2] <=> (split(/\t/,$b))[2] } @genes) {
  print $line;
}
