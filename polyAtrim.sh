#!/bin/bash
cmd="cat"
if [[ $1 =~ gz$ ]]; then cmd="zcat"; fi
l=10
if [ -n "$2" ]; then l=$2; fi
exe="$cmd $1 | fastq_to_fasta | sed 's/A\{$l\}.*$//'"
echo $exe; eval $exe
