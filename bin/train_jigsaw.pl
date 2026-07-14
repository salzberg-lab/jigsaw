#!/usr/bin/env perl

#################################################
## train_jigsaw.pl
##
## Purpose: Script to automate the training
##          of the statistical jigsaw
##
#################################################

use strict;
use Getopt::Long;
use File::Basename;

my $evidence_file;
my $ntree = 10;
my $treedir;
my $idir_list;
#my $adj_nt = 250;
## may want to make some of these command line args
## training evidence will be created in this file
## The purpose of this flag is to prepend to filenames the directory name to all
## file names.  Sometimes in each directory the files are prefixed with an identifier
## of the directory name
my $fasta;
my $tsplit;
my $addPrefix;
my $gridCmd;

my $HELP_INFO = q~
    train_jigsaw.pl [options]

    Input Options:

       -f <fasta file name>                 - the filename of the fasta formatted genomic sequence
       -e <evidence_file>                   - evidence file
       -d <output directory>                - directory to place decision trees in
       -l <input directories>               - the assumption is, a single directory contains
                                              all of the information for one genomic sequence
       -p <adjacent nt len>                 - specify how many nucleotides adjacent to gene model to
                                              include in training, defualt is 250.
       -n <number of decsion trees>         - number of decision trees to use (default is 10)
       -t <tree splits>                     - two types of decision trees can be used, when this option is set,
                                              the slower but potentially more accurate trees are used.
                                              The default is to use the method, which generates fast trees.
    Miscellaneous Options:
       -h, -help        - print this help message
       -V               - obtain program version
    ~;

my %Options;
my $err = GetOptions(\%Options,"e=s","n=s","d=s","l=s","p=s","f=s","x=s","t","help");
my $fasta = $Options{f};
my $evidence_file = $Options{e};
my $treedir = $Options{d};
my $idir_list = $Options{l};
my $adj_nt = 50;
$adj_nt = $Options{p} if($Options{p});
my $addPrefix = $Options{x};
my $tsplit = $Options{t};
my $helpval = $Options{help};
if($helpval) {
  print $HELP_INFO;
  exit(0);
}

$gridCmd="";
if ( defined $tsplit ) {
  $tsplit = "-a";
}

if( defined $addPrefix ) {
  $addPrefix = 0;
} else {
  $addPrefix = 1;
}

if ( ! defined $fasta ) {
    die("You must provide the filename of the fasta formatted genomic sequence expected in each directory  ");
}

if ( ! defined $idir_list ) {
    die("You must provide a file containing a list of each input directory using -l  ");
}

if ( ! defined $evidence_file ) {
    die("You must specify an evidence file using -e ");
}

if ( ! defined $treedir ) {
    #die("You must specify a directory to place the decision trees using -d");
  $treedir = ".";
}

my $tfname = $evidence_file;

my @asmbls;
open(FILE,$idir_list) || die("Cannot open $idir_list: $!");
while(my $line = <FILE>) {
  chomp($line);
  push( @asmbls, $line);
}
close(FILE);
my $evidence;
my $train_template;
open(FILE,$evidence_file) || die("Cannot open $evidence_file: $!");
while(my $line = <FILE>) {
  chomp($line);
  $evidence .= "$line:";
  if($line =~ /curation/ ) {
     my @tmp=split(/\s+/,$line);
     $train_template= $tmp[0];
  }
}
close(FILE);

if ( ! -d $treedir ) {
  system("mkdir $treedir");
}

my @vecIds = ("acc","don","code","start","stop","intron");
my @vlst;
for(my $iter = 0; $iter <= $#asmbls; $iter++) {
  my $dir = $asmbls[$iter];
  my $dir_id = basename($dir);
  my $prefix;
  if ( $addPrefix ) {
    $prefix = "$dir/${dir_id}_";
    if( ! (-e "${prefix}$fasta") ) {
    	$prefix = "$dir/$dir_id.";
    }
  } else {
    $prefix = "$dir/";
  }
  my $mytrain_file = $prefix . $train_template;
  #my $end5 = $adj_nt;
  #my $end3 = $adj_nt;
  my $vid = "ev_vec";
  $tfname = basename($tfname);
  writeTrainFile($tfname,$evidence,$prefix);
  my $vecfile = "$treedir/$dir_id.$vid";
  my $cmd = "jigsaw -f ${prefix}$fasta -e ${prefix}$tfname -t $vecfile -d $treedir -q $adj_nt";
  print "[$cmd]\n";
  system($cmd);
  remove_incomplete("$vecfile.start");
  remove_incomplete("$vecfile.stop");
  for (my $jiter = 0; $jiter <= $#vecIds; $jiter++) {
    my $vf = "$vecfile.$vecIds[$jiter]";
    if ( !$vlst[$jiter] ) {
      $vlst[$jiter] = $vf;
    } else {
      $vlst[$jiter] = $vlst[$jiter] . " " . $vf;
    }
  }
}

my @con_data;
my @dtFiles;
for(my $iter = 0; $iter <= $#vecIds; $iter++) {
  my $vecs = convert_vectors_to_triplet_vectors($vlst[$iter]);
  my $cntOcc = 1;
  if( $iter == 2) {
    $cntOcc = 0;
  } 
  $con_data[$iter] = wght($vecs,$cntOcc);
  my $dtFile = writeFilesAndMakeTrees($treedir,$con_data[$iter],$vecIds[$iter]);
  push(@dtFiles,$dtFile);
}

foreach my $dtFile (@dtFiles) {
  gtree($gridCmd,$dtFile,$ntree,$tsplit);
}

my $cmd = "rm -f $treedir/*.ev_vec.*";
print "$cmd\n";
system($cmd);

sub writeFilesAndMakeTrees {
  my ($treedir,$data,$id) = @_;
  my $dtfile = "$treedir/$id";
  my $cnfile = "$treedir/$id.con";
  open(DTFILE,">$dtfile") || die "unable to write to [$dtfile]\n";
  open(CNFILE,">$cnfile") || die "unable to write to [$cnfile]\n";
  my @vals = split(/\n/,$data);
  for(my $iter = 0; $iter <= $#vals; $iter++) {
    my @subvals = split(/ /,$vals[$iter]);
    my $class = 1;
    if( $subvals[$#subvals] >= 0.5 ) {
      $class = 2;
    }
    for(my $si = 0; $si < $#subvals; $si++) {
      print DTFILE "$subvals[$si] ";
    }
    print DTFILE "$class\n";
    print CNFILE "$vals[$iter]\n";
  }
  close(DTFILE);
  close(CNFILE);
  return $dtfile;
}

sub gtree {
  my ($gridcmd,$file,$num_trees,$split) = @_;
  for(my $cnt = 1; $cnt <= $num_trees; $cnt++) {
    my $seedcmd = "date +%s |";
    open(TPIPE,$seedcmd) || die "unable to open [$seedcmd]\n";
    my $seed = <TPIPE>;
    $seed =~ s/\n//;
    close(TPIPE);
    my $cmd="";
    if( $gridcmd ) {
      my $output = "$file.dt.$cnt";
      $cmd = "${gridcmd} -c \"mktree -t $file -s $seed -D $output $split\" -nowait";
    } else {
      $cmd = "mktree -t $file -s $seed $split";
    }
    print "$cmd\n";
    system($cmd);
    if ( !$gridcmd ) {
      if ( -e "$file.dt" ) {
	$cmd = "mv $file.dt $file.dt.$cnt";
	system($cmd);
	if ( -e "$file.dt.unpruned" ) {
	  $cmd = "rm -f $file.dt.unpruned";
	  system($cmd);
	}
      } else {
	print "some unanticipated problems running mktree\n";
	print "[$cmd] did not create: [$file.dt]\n";
      }
      sleep(5);
    }
  }
  # don't need the dt file anymore
  if(!$gridcmd) {
    my $cmd = "rm -f $file";
    system($cmd);
  }
}

sub loadGenes {
  my ($cdfile) = @_;
  my %ranges;
  #my $cdfile = "$asmblId/$asmblId.ceres.curation.cd";
  if( ! -e $cdfile || -z $cdfile) {
    print "could not open training file [$cdfile]\n";
    return %ranges;
  }
  open(FILE,$cdfile) || die "unable topen, $cdfile";
  while(my $line = <FILE>) {
    chomp($line);
    my @vals = split(/ /,$line);
    my $end5Idx = 1;
    my $end3Idx = 2;
    if( $#vals == 4 ) { 
      ($end5Idx,$end3Idx) = (2,3);
    }
    my $id = $vals[0];
    if( $vals[$end5Idx] > $vals[$end3Idx] ) {
	 ($vals[$end5Idx],$vals[$end3Idx]) = ($vals[$end3Idx],$vals[$end5Idx]);
    }
    if( $ranges{$id} ) {
      my @curr = split(/ /,$ranges{$id});
      if( $curr[0] > $vals[$end5Idx] ) {
	$curr[0] = $vals[$end5Idx];
      }
      if( $curr[1] < $vals[$end3Idx] ) {
	$curr[1] = $vals[$end3Idx];
      }
      $ranges{$id} = $curr[0] . " " . $curr[1];
      #print ":[$ranges{$id}]\n";
    } else {
      $ranges{$id} = $vals[$end5Idx] . " " . $vals[$end3Idx];
      #print "[$ranges{$id}]\n";
    }
  }
  close(FILE);
  return %ranges;
}

sub writeTrainFile {
  my ($tfname,$evlst,$prefix) = @_;
  my @evlst = split(/:/,$evlst);
  my $traindata = "$prefix$tfname";
  print "write: $traindata\n";
  open(TRAINDATA,">$traindata") || die "unable to write [$traindata]\n";
  for(my $iter = 0; $iter <= $#evlst; $iter++) {
    print TRAINDATA "$prefix$evlst[$iter]\n";
  }
  close(TRAINDATA);
}

sub convert_vectors_to_triplet_vectors {
  my ($lst) = @_;
  my @files = split(/ /,$lst);
  my %lookup;
  for(my $iter = 0; $iter <= $#files; $iter++) {
#    print "red file: $files[$iter]\n";
    compress($files[$iter],\%lookup);
  }
  return printLookup(\%lookup);
}

#####################################
## assume come in ascending order
#####################################

sub compress {
  my ($file,$lookup_ref) = @_;
  #print "open $file\n";
  my $cmd = "sort -k1n,1n $file |";
  print "$cmd\n";
  if( !open(FILE,$cmd) ) {
    print "warning unable to open [$file]\n";
    return;
  }
  my @lines = ("","","","","","");
  my $fhalf = ($#lines+1)/2;
  my $fail = 0;
  for(my $iter = 0; $iter <= $fhalf-1; $iter++) {
      #print "$iter, $fhalf\n";
      if($lines[$iter] = <FILE>) {
	chomp($lines[$iter]);
      } else {
	$fail = 1;
	last;
      }
      if($lines[$iter+$fhalf] = <FILE>) {
	chomp($lines[$iter+$fhalf]);
      } else {
	$fail = 1;
	last;
      }
  }
  if( $fail ) {
    close(FILE);
    return;
  }
  while($lines[5]) {
      my $tvals="";
      for(my $iter = 0; $iter <= $fhalf - 1; $iter++) {
	$tvals .= $lines[$iter];
	if( $iter < $fhalf-1 ) {
	  $tvals .= "#";
	}
      }
      my ($id1,$answer1,$vid1,$isValid1,$strnd1) = makeId($tvals);
      if( $id1  ) {
      my $tvals2="";
      for(my $iter = $fhalf; $iter <= $#lines; $iter++) {
	$tvals2 .= $lines[$iter];
	if( $iter < $#lines ) {
	  $tvals2 .= "#";
	}
      }
      my ($id2,$answer2,$vid2,$isValid2,$strnd2) = makeId($tvals2);
      if( $isValid1 ) {
	doStore($id1,$answer1,$vid1,$strnd1,$lookup_ref);
      }
      if( $isValid2 ) {
	doStore($id2,$answer2,$vid2,$strnd2,$lookup_ref);
      }
      } else {
	#print "try to continue: [$tvals]\n";
      }
      #3
      for(my $iter = 0; $iter <= $fhalf - 2; $iter++) {
	$lines[$iter] = $lines[$iter+1];
      }
      if( $lines[$fhalf-1] = <FILE> ) {
	chomp($lines[$fhalf-1]);
      } else {
	last;
      }
      for(my $iter = $fhalf; $iter <= $#lines - 1; $iter++) {
	$lines[$iter] = $lines[$iter+1];
      }
      if( $lines[$#lines] = <FILE> ) {
	chomp($lines[$#lines]);
      } else {
	last;
      }
    }
  close(FILE);
}

sub printLookup {
  my ($lookup_ref) = @_;
  my $vals;
  my $cnt = 0;
  my $tstop = 0;
  foreach my $id (keys %$lookup_ref) {
    if( !$id ) {
      next;
    }
    my ($nocode, $code);
    if($$lookup_ref{$id}{0}) {
      $nocode =$$lookup_ref{$id}{0};
    } else {
      $nocode =0;
    }
    if($$lookup_ref{$id}{1}) {
      $code =$$lookup_ref{$id}{1};
    } else {
      $code =0;
    }
	my @vals = split(/ /,$id);
	if(!$nocode && !$code) {
	    print "problems...";
	    print "[$id]\n";
	    exit(1);
	}
	#my $exp = $nocode+$code;
	my $cnt = 0;
	for(my $iter = 0; $iter <= $#vals; $iter++) {
	  $cnt += $vals[$iter];
	}
	if( !$cnt ) {
	    next;
	}
	#print "$id ";
        $vals .= "$id ";
	$vals .= printCode($nocode);
	#print " ";
	$vals .= " ";
	$vals .= printCode($code);
	$vals .= "\n";
	#print "\n";
	$cnt++;
    }
  return $vals;
}

sub printCode {
    my ($string) = @_;
    my @vals = split(/\+/,$string);
    my $max = 0;
    my $sum = 0;
    my $size = scalar(@vals);
    for( my $iter = 0; $iter < $size; $iter++) {
      if( $vals[$iter] > $max ) {
	$max = $vals[$iter];
      }
      $sum += $vals[$iter];
    }
    my $avg = 0;
    if( $size > 0 ) {
	$avg = ($sum/$size);
    }
    if( $size == 1 && $vals[0] == 0) {
	$size = 0;
    }
    return "$max $size $sum $avg";
}

sub makeId {
  my ($lns) = @_;
  #print "makeId input: [$lns]\n";
  my @lines = split(/\#/,$lns);
  my $id = "";
  my $vid = "";
  my $answer = -1;
  my $isValid = 1;
  my $middle = $#lines/2;
  if( $middle != 1 || $#lines != 2 ) {
    print "error check failed, $middle, [$lns], $#lines\n";
    exit(1);
  }
  my $line0 = $lines[0];
  my @vals0 = split(/ /,$line0);
  my $line1 = $lines[1];
  my @vals1 = split(/ /,$line1);
  my $line2 = $lines[2];
  my @vals2 = split(/ /,$line2);
  my $lGap = $vals1[0]-$vals0[1];
  my $rGap = $vals2[0]-$vals1[1];
  my ($mEnd5,$mEnd3) = ($vals1[0],$vals1[1]);
  if( $lGap < 1 ) {
    #print "error: $lns, [$lGap] \n";
    return ("","","","","");
    #exit(1);
  }
  if( $rGap < 1 ) {
    #print "error 2 $rGap < 1, input $lns\n";
    return ("","","","","");
    #exit(1);
  }
  my $strnd;
  for(my $iter = 0; $iter <= $#lines; $iter++) {
    my $line = $lines[$iter];
    my @vals = split(/ /,$line);
    if(0) { # is there a good reason to keep these out?
    if( $iter == $middle ) {
      my $doCont = 0;
      for(my $cnt =3; $cnt <= $#vals-1; $cnt++) {
	if( $vals[$cnt] ) {
	  $doCont = 1;
	  last;
	}
      }
      if( ! $doCont ) {
	$isValid = 0;
	last;
      }
    }
  }
    my ($v1,$v2);
    if( $iter == 0 && $lGap > 1 ) {
      $v1=$v2=$mEnd5-1;
    } elsif( $iter == 2 && $rGap > 1 ) {
      $v1=$v2=$mEnd3+1;
    } else {
      ($v1,$v2) = ($vals[0],$vals[1]);
    }
    if( !$vid ) {
      $vid  = $v1 . " " . $v2;
    } else {
      $vid  = $vid . " " . $v1 . " " . $v2;
    }
    my $bp = $v2-$v1 + 1;
    my $subid = "";
    for(my $jiter = 3; $jiter < $#vals; $jiter++) {
      if( ($lGap == 1 && $iter==0) || $iter==1 || ($iter==2 && $rGap==1) ) {
	$subid .= $vals[$jiter];
      } else {
	$subid .= "0"
      }
      if( $jiter < $#vals - 1 ) {
	$subid .= " ";
      }
    }
    $strnd = $vals[2];
    if( !$id ) {
      $id = $subid;
    } elsif( $vals[2] eq "+" ) {
      $id = $id . " " . $subid;
    } else {
      $id = $subid . " " . $id;
    }
    if( $iter == $middle ) {
      $answer = $vals[$#vals];
    }
  }
  if( (!$id || $answer == -1) && $isValid ) {
    print "fatal error, [$id], [$answer], [$lns]\n";
    exit(1);
  }
  return ($id,$answer,$vid,$isValid,$strnd);
}

sub doStore {
  my ($id,$answer,$vls,$strnd,$lookup_ref) = @_;
  #print "[$id]\n [$answer]\n[$vls]\n";
  my @vals = split(/ /,$vls);
  my $middle = ($#vals+1)/2 - 1;
  my $bp = 0;
  for(my $iter = 0; $iter < $#vals; $iter+=2) {
    $bp += $vals[$iter+1]-$vals[$iter]+1;
  }
  my $v1 = $vals[$middle]-$vals[$middle-1];
  #print "$vals[$middle] $vals[$middle-1] $vals[$middle+2] $vals[$middle+1]\n";
  if( $v1 < 1 ) {
      print "error @vals $middle\n";
      exit(1);
  } else {
    #$v1 = 1.0 - (1.0/$v1);
  }
  my $v2 = $vals[$middle+2]-$vals[$middle+1];
  if( $v2 < 1 ) {
      print "error 2: @vals, $middle\n";
      exit(1);
  }
  my $diststr;
  if( $strnd eq "+" ) {
    $diststr = $v1 . " " . $v2;
  } else {
    $diststr = $v2 . " " . $v1;
  }
  if( $$lookup_ref{$id}{$answer} ) {
    $$lookup_ref{$id}{$answer} .=  "+" . $bp;
    #print "pow: $$lookup_ref{$id}{$answer}\n";
  }  else {
    #print "$id $answer\n";
    $$lookup_ref{$id}{$answer} = $bp;
    #print "poww: $$lookup_ref{$id}{$answer}\n";
  }
}

sub wght {
  my ($data,$cntOcc) = @_;
  my @data_arr = split(/\n/,$data);
  my $noCode = 0;
  my $code = 0;
  my $numFeat;
  my $noCodeIdx;
  my $codeIdx;
  my $results;

  for(my $iter = 0; $iter <= $#data_arr; $iter++) {
    my $line = $data_arr[$iter];
    my @vals = split(/ /,$line);
    $numFeat = $#vals-7;
    if( $cntOcc ) {
      $noCodeIdx = $#vals-6;
    } else {
      $noCodeIdx = $#vals-5;
    }
    $codeIdx = $noCodeIdx+4;
    $noCode += $vals[$noCodeIdx];
    $code += $vals[$codeIdx];
  }

  my $max = 0;
  for(my $iter = 0; $iter <= $#data_arr; $iter++) {
    my $line = $data_arr[$iter];
    $line =~ s/\n//;
    my @vals = split(/ /,$line);
    my $val1 = 0;
    if( $noCode > 0 ) {
      $val1 = ($vals[$noCodeIdx]/$noCode);
    }
    my $val2 = 0;
    if( $code > 0 ) {
      $val2 = ($vals[$codeIdx]/$code);
    }
    my $tval = 0;
    if( $val1+$val2 > 0 ) {
      $tval = -log($val1+$val2);
    }
    if($tval > $max) {
      $max = $tval;
    }
  }

  if( $max <= 0) {
    print "training error - probably a directory structure/file format isue.\n";
    exit(1);
  }

  for(my $iter = 0; $iter <= $#data_arr; $iter++) {
    my $line = $data_arr[$iter];
    $line =~ s/\n//;
    my @vals = split(/ /,$line);
    my $val1 = 0;
    if( $noCode > 0 ) {
      $val1 = ($vals[$noCodeIdx]/$noCode);
    }
    my $val2 = 0;
    if( $code > 0 ) {
      $val2 = ($vals[$codeIdx]/$code);
    }
    my $tnc  = $vals[$noCodeIdx];
    my $tc = $vals[$codeIdx];
    my $prob = 0;
    if( $tnc+$tc > 0 ) {
      $prob = $tc/($tnc+$tc);
    }
    my $numPred = 0;
    for(my $iter = 0; $iter < $numFeat; $iter++) {
      if($vals[$iter] > 0) {
	$numPred++;
      }
    }
    for(my $iter = 0; $iter < $numFeat; $iter++) {
      $results .= "$vals[$iter] ";
    }
    my $tval = 0;
    if( $val1+$val2 > 0 ) {
      $tval = -log($val1+$val2);
    }
    if( $max > 0 ) {
      $tval = ($tval/$max);
    }
    $tval = 1-$tval;
    $tval += 0.01;
    $results .= "$prob\n";
      #print "$prob $tval ";
      #print "\n";
  }
  return $results;
}

sub remove_incomplete {
  my ($file) = @_;
  open(FILE,$file) || die "unable to open [$file]\n";
  my $pass = 0;
  while(my $line = <FILE>) {
    chomp($line);
    my @vals = split(/ /,$line);
    if( $vals[$#vals] ) {
      $pass = 1;
      last;
    }
  }
  close(FILE);
  if( !$pass ) {
    print "looks like incomplete gene removing file....$file\n";
    system("rm -f $file");
  }
}
