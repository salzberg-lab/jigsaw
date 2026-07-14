#!/usr/bin/perl -w

use strict;

my $root = $ARGV[0] || die "need the root directory\n";
my $ansFile=$ARGV[1] || die "must specify file containing gene coords\n";
my $tfid = $ARGV[2] || die "must specify the filename id\n";
my $odir = $ARGV[3] || die "must specify an output directory\n";
my $aid = $ARGV[4] || "";
if($aid) {
  $aid.= ".";
}

my $filt="";
my $pfid=$tfid;
if($tfid =~ / /) {
  ($pfid,$filt)=split(/ /,$tfid);
}

my $SYSRUN=1;

my $mgFile = "$odir/$aid$pfid.mg";
my $run1 = "miss.pl $root $ansFile 1 $pfid $filt > $mgFile";
print "$run1\n";
system($run1) if($SYSRUN);
my $mgR = lineCnt($mgFile);

my $meFile = "$odir/$aid$pfid.me";
my $run2 = "miss.pl $root $ansFile 0 $pfid $filt > $meFile";
print "$run2\n";
system($run2) if($SYSRUN);
my $meR = lineCnt($meFile);

my $weFile = "$odir/$aid$pfid.we";
my $run3 = "we.pl $root $ansFile $pfid $filt > $weFile";
print "$run3\n";
system($run3) if($SYSRUN);
my $weR = lineCnt($weFile);

my $ntfile="$odir/$aid$pfid.spec";
my $run4 = "spec.pl $root $ansFile $pfid $filt > $ntfile";
print "$run4\n";
system($run4) if($SYSRUN);
my $ofile="$odir/$aid$pfid.report";
open(FILE,">$ofile") || die "unable to write [$ofile]\n";
#print FILE "$ntgR $npgR $ncgR $sng $mgR $texR $cexR $cep $meR $weR $snR $spg $avgt $accRatio\n";
print FILE "whole $mgR $meR $weR\n";
open(FILE2,"$ntfile") || die "unable to open [$ntfile]\n";
while(my $line = <FILE2>) {
  if( $line =~ /^raw/ || $line =~ /^sng:/ ) {
    print FILE $line;
  }
}
close(FILE2);
close(FILE);

sub lineCnt {
  my ($file) = @_;
  my $cmd = "wc -l $file | ";
  open(PIP,$cmd) || die "unable to run [$cmd]\n";
  my $ln=0;
  my $line=<PIP>;
  if($line=~ /\s*(\d+)/ ) {
    $ln = $1;
    #print "lc: $ln\n";
  }
  close(PIP);
  return $ln;
}

sub parseEvalFile {
  my ($file)=@_;
  open(FILE,$file) || die "unable to open [$file]\n";
  while(my $line=<FILE>) {
    chomp($line);
    if( $line=~ /Summary/ ) {
      my @vals=split(/\s+/,$line);
      my ($ntg,$npg,$ncg,$sn,$tex,$cex) = ($vals[2],$vals[3],$vals[4],$vals[15],$vals[20],$vals[22]);
      return ($ntg,$npg,$ncg,$sn,$tex,$cex);
    }
  }
  close(FILE);
  return "";
}
