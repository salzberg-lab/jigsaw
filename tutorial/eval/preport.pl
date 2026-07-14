#!/usr/bin/perl -w

use strict;

my $id;
printf("%38s %5s %5s %5s %5s %6s %5s %6s %6s %6s %6s %5s %5s %5s\n",
       "PROG","SNG", "SPG","SNE","SPE","SNI","SPI","SNN","SPN","PGC","TGC","MG","ME","WE");
my ($mg,$me,$we,$te,$ti);
my @results;
my ($pg,$tg);
while(my $line=<STDIN>) {
  chomp($line);
  if( $line =~ /==>(.*)<==/ ) {
    $id = $1;
    $id =~ s/\s+$//g;
    next;
  } elsif( $line =~ /^::::::/ ) {
    $id = <STDIN>;
    chomp($id);
    <STDIN>;
    next;
  } elsif( $line =~ /^raw / ) {
    $line =~ s/^raw //g;
    my @vals = split(/ /,$line);
    ($pg,$tg,$te,$ti) = ($vals[0],$vals[8],$vals[9],$vals[10]);
    next;
  } elsif( $line =~ /^whole/ ) {
    $line =~ s/^whole //g;
    ($mg,$me,$we)=split(/\s+/,$line);
  } elsif( $line =~ /^sng/ ) {
    my @vals=split(/\s+/,$line);
    #print "$id $line\n";
    if( !$id ) {
      $id="prog";
    }
    my $gc=($vals[1]+$vals[3])/2;
    my $ec=($vals[5]+$vals[7])/2;
    my $ic=($vals[9]+$vals[11])/2;
    my $nc=($vals[13]+$vals[15])/2;
    my $fl=sprintf("%38s  %2.3f %2.3f  %2.3f %2.3f  %2.3f %2.3f  %2.3f %2.3f  %5d %5d  %2.3f %2.3f %2.3f",
		   $id, $vals[1],$vals[3],$vals[5],$vals[7],$vals[9],$vals[11],$vals[13],$vals[15],$pg,$tg,$mg/$tg,$me/$te,$we/$ti);
    push (@results,$fl);
  }
}

my @sres = sort { (split(/\s+/,$a))[2] <=> (split(/\s+/,$b))[2] } @results;
foreach my $val (reverse @sres) {
  print "$val\n";
}

print "\nTable Key\n\n";
print "PROG = program name\n";
print "SNG = CG/TG, ";
print "SPG = CG/PG, ";
print "SNE = CE/TE, ";
print "SPE = CE/PE, ";
print "SNI = CI/TI, ";
print "SPI = CI/PI, ";
print "SNN = CN/TN, ";
print "SPN = CN/PN\n";
print "PGC = Predicted Gene Count\n";
print "TGC = Known True Gene Count\n";
print "MG = Completely Missed Genes\n";
print "ME = Completely Missed Exons\n";
print "WE = Whole Wrongly Inserted Exons\n";
print "\n";
print "PG = Predicted Gene Count\n";
print "PE = Predicted Exon Count\n";
print "PI = Predicted Intron Count\n";
print "PN = Predicted Protein Coding Nucleotides\n";
print "\n";
print "CG = Correctly Predicted Genes\n";
print "CE = Correctly Predicted Exons\n";
print "CI = Correctly Predicted Introns\n";
print "CN = Correctly Predicted Protein Coding Nucleotides\n";
print "\n";
print "TG = Known True Genes\n";
print "TE = Known True Exons\n";
print "TI = Known True Introns\n";
print "TN = Known True Protein Coding Nucleotides\n\n";
