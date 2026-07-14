#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;

my $HELP_INFO = q~
    run_combiner.pl [options]

    Input Options:
       -f <fasta file> - name of fasta formatted genomic sequence
       -o <combiner output file name> - name of file for gene predictions
       -e <evidence_file>    - template of evidence file, should be same one used for training
       -d <output directory>    - directory to read decision trees from
       -l <input directories>    - the assumption is, a single directory contains
                                   all of the information for one genomic sequence
                                   this is a file containing a list of those directories (1 per line)
	    -n <intron penalty> - specify maximum intron length and penalty for exceeding the length
	    -i <minimum intron length> -  specify minimum intron length
       -c  - do not overwrite existing output file (default is to overwrite)
       -lin - run linear combiner
       -condor 							- run on condor machine 
    Miscellaneous Options:
       -h, -help        - print this help message
       -V               - obtain program version
    ~;

my %Options;
my $err = GetOptions(\%Options,"e=s","o=s","d=s","l=s","f=s","x","lin","n=s","i=s","c","help");
my $output = $Options{o};
my $myfasta = $Options{f};
my $evidence_file = $Options{e};
my $treedir = $Options{d};
my $idir_list = $Options{l};
my $addPrefix = $Options{x};
my $useLinComb = $Options{lin};
my $noOverWrite = $Options{c};
my $penStr = $Options{n};
my $inLen = $Options{i};
my $helpval = $Options{help};

my $localOnly=1;
if ( defined $addPrefix ) {
  $addPrefix = 0;
} else {
  $addPrefix = 1;
}

if ( defined $helpval ) {
  print $HELP_INFO;
  exit(0);
}

if ( ! defined $myfasta ) {
    die ("You must specify a fasta file using -f ");
}

if ( ! defined $output ) {
    die ("You must specify an file name for Combiner output using -o ");
}

if ( ! defined $idir_list ) {
    die("You must specify file listing directory names using -l ");
}

if ( ! defined $evidence_file ) {
    die("You must specify an evidence file using -e ");
}

if ( ! defined $treedir && ! defined $useLinComb ) {
    die("You must specify a directory to place the decision trees using -d");
}

my $penArg="";
if ( defined $penStr ) {
   $penArg="-n \"$penStr\"";
}

my $linArg="";
if ( defined $useLinComb ) {
   $linArg = "-l";
}

my $ilArg="";

my $evidence;
open(FILE,$evidence_file) || die("Cannot open $evidence_file: $!");
while (my $line = <FILE>) {
  chomp($line);
  ## ignore the "answers line, not used when actually running combiner
  if ( $line !~ / curation/ ) {
    $evidence .= "$line:";
  }
}
close(FILE);

my @testOn;
open(FILE,$idir_list) || die "unable to open [$idir_list]\n";
while (my $line = <FILE>) {
  chomp($line);
  push(@testOn,$line);
}
close(FILE);
my $combid = basename($evidence_file);
for (my $cnt = 0; $cnt <= $#testOn; $cnt++) {
  my $dir = $testOn[$cnt];
  my $dname = basename($dir);
  my $prefix;
  if ( $addPrefix ) {
    $prefix = "$dir/$dname.";
  } else {
    $prefix = "$dir/";
  } 
  my $treearg="";
  $treearg = "-d $treedir" if(defined $treedir);
  my $combd = "${prefix}${combid}";
  my $succ = writeEvidenceFile($combd,$evidence,$prefix);
  if ($succ) {
    my $combout = "${prefix}$output";
    my $cmd;
    if ( !$localOnly ) {
      $cmd = "condor_run \"jigsaw -f ${prefix}$myfasta $treearg $linArg -e $combd -m $combout $penArg\" >& /dev/null&";
    } else {
      $cmd = "jigsaw -f ${prefix}$myfasta $treearg $linArg -e $combd -m $combout $penArg";
    }
    print "[$cmd]\n";
    if ( !defined($noOverWrite) || ! -e $combout || -z $combout ) {
      system($cmd);
    }
  }
}

sub writeEvidenceFile {
  my ($tfname,$evlst,$prefix) = @_;
  my @evlst = split(/:/,$evlst);
  my $evdata = "$tfname";
  my $status = 0;
  if( open(EVDATA,">$evdata") ) {
     for(my $iter = 0; $iter <= $#evlst; $iter++) {
       print EVDATA "$prefix$evlst[$iter]\n";
     }
     close(EVDATA);
     $status = 1;
  } else {
     die "unable to write [$evdata]\n";
  }
  return $status;
}
