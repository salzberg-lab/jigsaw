#!/usr/bin/perl

#$infile = $ARGV[0];

&Parse_GlimmerM_out();

sub Parse_GlimmerM_out {
    local($infile)=@_;
    local(@x,@a,$line,$i,$j);
    $DEBUG = 1;
    #open(G,"$infile") || die "Can't open Glimmer output $out\n";

    $cnts{1} = 0;
    $cnts{2} = 0;
    while($line = <STDIN>) {
      $line =~ s/\n//;
      if(!$line) {
	next;
      }
      if ($line =~ /^
            \s*(\d+)
	    \s+(\d+)  ## first field:        \1
                           ##   gene_number.exon_number 
            \s+([+-])        ## Strand (+, input strand; - opposite) \3
	  \s+(\w+)
            \s+(\d+)       ## End5  \3
            \s+(\d+)       ## End3  \4
            \s+\d+         ## Length
            /x) {
	($model, $num, $orient, $type , $end5, $end3) = ($1, $2, $3, $4, $5, $6);
	my $score = 1.0;
	my $frame=1.0;
	print "seqid\tglm\t$type\t$end5\t$end3\t$score\t$orient\t$frame\t$model\n";
      } else {
      }
    }
    close(G);
}
