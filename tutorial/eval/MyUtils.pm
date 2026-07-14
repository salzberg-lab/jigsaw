package MyUtils;

use Exporter;
@ISA=('Exporter');
@EXPORT=qw(&readFasta &writeFasta &compRel &revComplement &compHist &printHist &shuffle);

sub shuffle {
  my $array = shift;
  my $i;
  for ($i = @$array; --$i; ) {
    my $j = int rand ($i+1);
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }
}

sub revComplement {
  my ($seq) = @_;
  $seq =~ tr/ACGTUuacgtyrkm/TGCAAatgcarymk/;
  $seq = reverse($seq);
  return $seq;
}

sub compRel {
  my ($crd,$strnd,$beg,$end) = @_;
  my $nc="";
  if ( $strnd eq "-") {
    $nc = ($end-$crd)+1;
  } elsif ($strnd eq "+") {
    $nc = ($crd-$beg)+1;
  }
  return $nc;
}

sub readFasta {
  my ($file)=@_;
  my %seqs;
  open(FILE,$file) || die "failed to read [$file]\n";
  my $header=<FILE>;
  if( !$header ) {
    return %seqs;
  }
  chomp($header);
  $header =~ s/^>//g;
  $header =~ s/\s+$//g;
  for (;;) {
    my $line=<FILE>;
    chomp($line) unless !$line;
    #print "line [$line]\n";
    if ( !$line || $line =~ />/ ) {
      #print "hd: $header\n";
      $seqs{$header} = $seq;
      $seq="";
      $header=$line;
      $header =~ s/^>//g unless !$line;
      $header =~ s/\s+$//g unless !$line;
    }
    if ( !$line ) {
      last;
    }
    if ( $line !~ />/ ) {
      $seq .= $line;
    }
  }
  close(FILE);
  return %seqs;
}


sub writeFasta {
  my ($str,$header,$fh) = @_;
  print $fh ">$header\n";
  #print ">$header\n";
  for(my $iter = 0; $iter < length($$str); $iter += 80) {
    my $mstr = substr($$str,$iter,80);
    print $fh "$mstr\n";
  }
}


sub compHistDeleteMe {
  my ($aRef,$dBins,$min,$max) = @_;

  my @vals = @$aRef;

  if( ! defined($min) || !defined($max) ) {
    my @svals=sort {$a<=>$b} @vals;
    $min=$svals[0];
    $max=$svals[$#svals];
  }
  my $range = ($max-$min)+1;
  if ( $dBins > $range ) {
    $range=$dBins;
  }
  my $modif=$range/$dBins;

  print "$min $max $range $modif $dBins\n";

  my %hist;

  my $tot=0;
  foreach my $val (@vals) {
    my $bin = int(($val-$min)/$modif);
    #print "$val, $bin\n";
    $hist{$bin}++;
    $tot++;
  }

  my @skeys = sort { $a <=> $b } keys %hist;
  foreach my $key (@skeys) {
    my $pcnt = $hist{$key}/$tot;
    print "$key $pcnt\n";
  }
}


sub compHist {
  my ($hRef,$min,$max,$nBins) = @_;
  if( !$min || !$max ) {
    my @svals=sort {$a<=>$b} keys %$hRef;
    $min=$svals[0];
    $max=$svals[$#svals];
  }
  if( !$min || !$max ) {
    print "what happend [$min][$max][$nBins]\n";
    exit(1);
  }
  my %hist;
  my $tot=0;
  my $bSize = (($max-$min)+1)/$nBins;
  print "$min $max $nBins $bSize\n";
  foreach my $val (keys %$hRef) {
    my $bin=int(($val-$min)/$bSize);
    my $cnt=$hRef->{$val};
    $hist{$bin}+=$cnt;
    #print "$val $min\n";
    #print "bin: $bin $hist{$bin}\n";
    $tot+=$cnt;
  }
  my @svals = sort { $a <=> $b } keys %hist;
  my @result;
  for(my $iter=0; $iter < $nBins; $iter++) {
    $result[$iter]=0;
  }
  foreach my $bin (@svals) {
    $result[$bin] = $hist{$bin}/$tot;
    #$result[$bin] = $hist{$bin};
    #print "$bin $result[$bin]\n";
  }
  return (\@result,$min,$bSize);
}

sub printHist {
  my ($fh,$aRef,$min,$step) = @_;
  my @arr = @$aRef;
  my $score = $min;
  for(my $iter=0; $iter <= $#arr; $iter++) {
    if( $arr[$iter] > 0 ) {
      if ( !$fh ) {
	print "$score $arr[$iter]\n";
      } else {
	#print "$score $arr[$iter]\n";
	print $fh "$score $arr[$iter]\n";
      }
    }
    $score += $step;
  }
  if( !$fh ) {
	print "\n\n";
  } else {
     print $fh "\n\n";
  }
} 

