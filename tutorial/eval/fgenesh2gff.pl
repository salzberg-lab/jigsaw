#!/usr/bin/perl -w

use strict;

my @genes;
while(my $line= <STDIN>) {
  chomp($line);
  #print "$line\n";
  if( $line =~ /\s+(\d+) ## geneid
                \s+([-+]) ## strnd
                \s+(\d+) ## exon number
                \s+(\w{4}) ## exon type
                \s+(\d+) ## end5
                \s+\-
                \s+(\d+) ## end3
                \s+?(\d+\.\d+) ## score
      .*/x ) {
    my ($gid,$strnd,$exNum,$exonId,$end5,$end3,$score) = ($1,$2,$3,$4,$5,$6,$7);
    my $frame=0;
    push (@genes, "seqid\tfgenesh\t$exonId\t$end5\t$end3\t$score\t$strnd\t$frame\t$gid\n");
  }
}

#@genes = map { $_- }
         #map { $a->[1] <=> $b->[1] }
         #map { [$_, (split/\t/)[3] ] }
         #@genes;

foreach my $line (sort { (split(/\t/,$a))[3] <=> (split(/\t/,$b))[3] } @genes) {
  print "$line";
}

