#!/bin/sh -xvf

root=../
evalf=$root/prog_eval/ans.txt
outdir=$root/prog_eval
file=$1
report.pl $root $evalf "$file" $outdir
