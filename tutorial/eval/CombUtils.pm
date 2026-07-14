package CombUtils;

use Exporter;

@ISA=('Exporter');
@EXPORT=qw(&loadTrueGenes &loadGff &loadGff2 &loadTrueGenes2 &loadDefault &loadGffHash &checkOrf);


sub loadTrueGenes {
  my ($g) = @_;
  open(F,$g) || die "unable to open $g for reading\n";

  my %truegenes;
  my %seen;
  my $strnd;
  my $lastModel = 0;
  my ($idx1,$idx2) = (2,3);
  my $modelIdx = 1;
  while(<F>) {
    chomp;
    my @a=split;
    my $geneId = $a[0];
    if ( $a[$idx1] > $a[$idx2]) {
      ($a[$idx1],$a[$idx2]) = ($a[$idx2],$a[$idx1]);
      $strnd = "-";
    } else {
      $strnd = "+";
    }
    if (!$seen{$geneId}) {
      $seen{$geneId}++;
      $truegenes{$geneId}=$strnd." ".$a[$idx1]." ".$a[$idx2];
    } else {
      if ( !($a[$modelIdx] eq $lastModel) ) {
	$truegenes{$geneId}.=":".$strnd;
      }
      $truegenes{$geneId}.=" ".$a[$idx1]." ".$a[$idx2];
    }
    $lastModel = $a[$modelIdx];
  }
  close(F);
  return %truegenes;
}

sub loadGff {
  my ($file,$filt) =  @_;
  if ($filt) {
    my $cmd = "cat $file | $filt |";
    #print "cmd: [$cmd]\n";
    open(F,$cmd) || die "unable to run [$cmd]\n";
  } elsif ( !open(F,$file) ) {
    print "unable to open $file for reading\n";
    return "";
  }
  my $genes="";
  my $lastModel;
  my $strndIdx = 6;
  #my $modelIdx = 8;
  my ($idx1,$idx2) = (3,4);
	my $tcnt=0;
  while (my $line=<F>) {
    chomp($line);
    if ($line =~ /\#/ ) {
      $lastModel="";
      next;			# skip coment line
    }
    $line =~ s/		/	/g;
    my @vals = split(/\s+/,$line);
	 my $modelIdx = $#vals;
	#foreach my $v (@vals) {
	#print "bs: $v\n";
		#if( !defined($v) ) {
			#print "failure: [$v]\n";
			#exit(1);
		#}
   #}
	#exit(1);
    my $strnd = $vals[$strndIdx];
	 #print "[$vals[$modelIdx]]\n";
    if (!defined($lastModel) ) {
		$tcnt++;
      $genes = $strnd;
		if($tcnt > 1 ){
       print "crap[$line]\n";
			exit(1);
		}
    }
    if ( defined($lastModel) &&  !($vals[$modelIdx] eq $lastModel) ) {
      $genes .= ":".$strnd;
    }
    #print "[@vals]\n";
    $genes .= " " . $vals[$idx1] . " " . $vals[$idx2];
    $lastModel = $vals[$modelIdx];
  }
  close(F);
  return $genes;
}


sub loadGff2 {
  my ($file,$filt) =  @_;
  if ($filt) {
    my $cmd = "cat $file | $filt |";
    #print "cmd: [$cmd]\n";
    open(F,$cmd) || die "unable to run [$cmd]\n";
  } elsif ( !open(F,$file) ) {
    print "unable to open $file for reading\n";
    return "";
  }
  my %genes=();
  my $lastModel;
  my $strndIdx = 6;
  my $modelIdx = 8;
  my ($idx1,$idx2) = (3,4);
  while (my $line=<F>) {
    chomp($line);
    if ($line =~ /\#/ ) {
      $lastModel="";
      next;			# skip coment line
    }
    my @vals = split(/\t/,$line);
    my $seqname = $vals[0];
    my $strnd = $vals[$strndIdx];
    if (!defined($lastModel) ) {
      $genes{$seqname} = $strnd;
    }
    if ( defined($lastModel) &&  !($vals[$modelIdx] eq $lastModel) ) {
      $genes{$seqname} .= ":".$strnd;
    }
    #print "[@vals]\n";
    $genes{$seqname} .= " " . $vals[$idx1] . " " . $vals[$idx2];
    $lastModel = $vals[$modelIdx];
  }
  close(F);
  return %genes;
}

sub loadGffHash {
  my ($file,$filt) =  @_;
  if ($filt) {
    my $cmd = "cat $file | $filt |";
    #print "cmd: [$cmd]\n";
    open(F,$cmd) || die "unable to run [$cmd]\n";
  } elsif ( !open(F,$file) ) {
    print "unable to open $file for reading\n";
    return "";
  }
  my %genes=();
  my $strndIdx = 6;
  my $modelIdx = 8;
  my ($idx1,$idx2) = (3,4);
  while (my $line=<F>) {
    chomp($line);
    if ($line =~ /\#/ ) {
      next;       # skip coment line
    }
    my @vals = split(/\t/,$line);
    my $strnd = $vals[$strndIdx];
	 my $geneId = $vals[$#geneId];
    if (!defined($genes{$geneId}) ) {
      $genes{$geneId} = $strnd;
    } else {
      $genes{$geneId} .= ":$strnd";
    }
    $genes{$geneId} .= " " . $vals[$idx1] . " " . $vals[$idx2];
  }
  close(F);
  return %genes;
}


sub loadTrueGenes2 {
  my ($g) = @_;
  open(F,$g) || die "unable to open $g for reading\n";

  my %truegenes;
  my %seen;
  my $strnd;
  my $lastModel = 0;
  my ($idx1,$idx2) = (2,3);
  my $modelIdx = 1;
  while(<F>) {
    chomp;
    my @a=split;
    my $seqId = $a[0];
    if ( $a[$idx1] > $a[$idx2]) {
      ($a[$idx1],$a[$idx2]) = ($a[$idx2],$a[$idx1]);
      $strnd = "-";
    } else {
      $strnd = "+";
    }
    my $geneId = $a[$modelIdx];
    if (!$truegenes{$seqId}{$geneId}) {
      $truegenes{$seqId}{$geneId}=$strnd." ".$a[$idx1]." ".$a[$idx2];
    } else {
      $truegenes{$seqId}{$geneId}.=" ".$a[$idx1]." ".$a[$idx2];
    }
  }
  close(F);
  return %truegenes;
}

sub loadDefault {
  my ($g) = @_;
  open(F,$g) || die "unable to open $g for reading\n";

  my $truegenes="";
  my $strnd;
  my ($idx1,$idx2) = (1,2);
  my $modelIdx = 0;
  my $lastGeneId="";
  while(<F>) {
    chomp;
    my @a=split;
    my $geneId = $a[$modelIdx];
    if ( $a[$idx1] > $a[$idx2]) {
      ($a[$idx1],$a[$idx2]) = ($a[$idx2],$a[$idx1]);
      $strnd = "-";
    } else {
      $strnd = "+";
    }
    if( !($geneId eq $lastGeneId) && $lastGeneId ) {
      $truegenes.= ":". $strnd." ".$a[$idx1]." ".$a[$idx2];
    } elsif( !$lastGeneId) {
      $truegenes=$strnd." ".$a[$idx1]." ".$a[$idx2];
    } else {
      $truegenes.=" ".$a[$idx1]." ".$a[$idx2];
    } 
    $lastGeneId =$geneId;
  }
  close(F);
  return $truegenes;
}

sub checkOrf {
  my ($exonStr,$seqRef,$checkOrf,$checkStartStop) = @_;
  my $seq = $$seqRef;
  my @exons = split(/ /,$exonStr);
  my $strnd = $exons[0];

  my $cnt = 0;
  my $errors="";
  for(my $iter = 1; $iter <= $#exons - 2; $iter += 2) {
    my $end3 = $exons[$iter+1];
    my $splice = lc(substr($seq,$end3+1-1,2));
    if(($strnd eq "-" && $splice eq "ct") || ($strnd eq "+" && ($splice eq "gt" || $splice eq "gc"))) {
    } else {
      my $type = "don";
      if( $strnd eq "-" ) {
	$type="acc";
      }
      $errors .= "end3: $end3, $type $splice bad\n";
    }
  }
  for(my $iter = 3; $iter <= $#exons - 1; $iter += 2) {
    my $end5 = $exons[$iter];
    my $splice = lc(substr($seq,$end5-2-1,2));
    if(($strnd eq "-" && ($splice eq "ac" || $splice eq "gc")) || ($strnd eq "+" && $splice eq "ag")) {
    } else {
      my $type = "acc";
      if( $strnd eq "-" ) {
	$type="don";
      }
      $errors .= "end5 $end5, $type $splice bad\n";
    }
  }
  if( $checkStartStop ) {
    my $end5 = $exons[1];
    my $str5 = lc(substr($seq,$end5-1,3));
    my $str51 = ""; #lc(substr($seq,$end5-4,3));
    if( ($strnd eq "+" && $str5 eq "atg") ||
	($strnd eq "-" && ($str5 eq "tca" || $str5 eq "cta" || $str5 eq "tta" || $str51 eq "tca" || $str51 eq "cta" || $str51 eq "tta")) ) {
    } else {
      my $type = "start";
      if( $strnd eq "-" ) {
	$type="stop";
      }
      $errors .= "end5 $end5, $type $str5 bad\n";
    }
    my $end3 = $exons[$#exons];
    my $str3 = lc(substr($seq,$end3-3,3));
    my $str31 = ""; #lc(substr($seq,$end3,3));
    if( ($strnd eq "+" && ($str3 eq "tga"||$str3 eq "tag"||$str3 eq "taa" || $str31 eq "tga"||$str31 eq "tag"||$str31 eq "taa"))
	|| ($strnd eq "-" && $str3 eq "cat") ) {
    } else {
      my $type = "stop";
      if( $strnd eq "-" ) {
	$type="start";
      }
      $errors .= "end3 $end3, $type $str3 bad\n";
    }
  }
  if ($checkOrf) {
    my $theCoding = "";
    my $sum = 0;
    for (my $iter = 1; $iter <= $#exons - 1; $iter += 2) {
      my $end5 = $exons[$iter];
      my $end3 = $exons[$iter+1];
      my $mlen = ($end3-$end5)+1;
      $theCoding .= substr($seq,$end5-1,($end3-$end5)+1);
      $sum += $mlen;
    }
    my $len = length($theCoding);
    if ( $len != $sum) {
      print "big problem here: $exonStr, len:[$len] sum:[$sum]\n";
      exit(1);
    }
    if ( $len % 3 != 0 ) {
      $errors .= "not multiple of 3: $len\n";
    }
    for (my $iter = 0; $iter < $len; $iter+=3) {
      my $codon = lc(substr($theCoding,$iter,3));
      if (isStop($codon,$strnd) && (($iter < $len-3 && $strnd eq "+") || ($iter > 2 && $strnd eq "-"))) {
	$errors .= "inframe stop $iter, $codon, ($len)\n";
      }
    }
  }
  return $errors;
}

sub isStop {
  my ($codon,$strnd) = @_;
  my $res=0;
  if( ($strnd eq "+" && ($codon eq "tga"||$codon eq "tag"||$codon eq "taa")) ||
      ($strnd eq "-" && ($codon eq "tca" || $codon eq "cta" || $codon eq "tta")) ) {
    $res=1;
  }
  return $res;
}
