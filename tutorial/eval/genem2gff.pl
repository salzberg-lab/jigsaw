#!/usr/bin/perl -w

use strict;

my $score = 1.0;
my @genes;

while(my $line = <STDIN>) {
  chomp($line);
  $line =~ s/^\s+//g;
  my @vals = split(/\s+/,$line);
  my $score = 1.0;
  my ($id,$type,$strnd,$end5,$end3) = ($vals[0],$vals[3],$vals[2],$vals[4],$vals[5]);
  if( $#vals >=5 && $vals[3] ) {
    #print "type? $type\n";
    #print "[$vals[0]] [$vals[1]] [$vals[2]]\n";
    if($type eq "Initial" ||
       $type eq "Internal" ||
       $type eq "Terminal" ||
       $type eq "Single" ) {
    #exit(1);
    push (@genes, "seqid\tgenemarkhmm\t$type\t$end5\t$end3\t$score\t$strnd\t.\t$id\n");
  }
  }
}

foreach my $line (sort { (split(/\t/,$a))[3] <=> (split(/\t/,$b))[3] } @genes) {
  print "$line";
}
