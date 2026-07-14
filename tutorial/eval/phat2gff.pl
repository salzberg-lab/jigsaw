#!/usr/bin/perl -w

use strict;

my $score = 1.0;
my @genes;

while(my $line = <STDIN>) {
  chomp($line);
  $line =~ s/^\s+//g;
  my @vals = split(/\s+/,$line);
  #my $score = 1.0;
  my ($id,$type,$strnd,$end5,$end3,$frame,$score) = ($vals[0],$vals[3],$vals[2],$vals[4],$vals[5],$vals[7],$vals[9]);
  if( defined($type) ) {
    #print "type? $type\n";
    #print "[$vals[0]] [$vals[1]] [$vals[2]]\n";
    if($type eq "init" ||
       $type eq "term" ||
       $type eq "intl" ||
       $type eq "sing" ) {
    #exit(1);
    push (@genes, "seqid\tgenemarkhmm\t$type\t$end5\t$end3\t$score\t$strnd\t$frame\t$id\n");
  }
  }
}

foreach my $line (sort { (split(/\t/,$a))[3] <=> (split(/\t/,$b))[3] } @genes) {
  print "$line";
}
