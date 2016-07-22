#! /usr/bin/env bash
set -e
if [ "$#" -lt 1 ]; then
	echo "USAGE: trimmomagic FASTQ [[FASTQ OR .] [[QUALITY] [[WINDOW] [[LENGTH] [CPU]]]]]" >&2
	echo "EXAMPLE: trimmomagic in.fq.gz . 25 3 6 8" >&2
	echo "EXAMPLE: trimmomagic in_1.fq in_2.fq 15 4 10 4" >&2
	exit 1
fi
in="$1"
first=$(echo "$1" | sed 's/.gz//;s/.fq//;s/.fastq//')
out="$first.qt.fq.gz"
mode='SE'
if [ "$#" -ge 2 ] && [ $2 != '.' ]; then
	in="$1 $2"
	second=$(echo "$2" | sed 's/.gz//;s/.fq//;s/.fastq//')
	out="$first.qt.fq.gz $first.nopair.fq.gz $second.qt.fq.gz $second.nopair.fq.gz"
	mode='PE'
fi
quality='25' && [ "$#" -lt 3 ] || quality=$3
window='3' && [ "$#" -lt 4 ] || window=$4
length='6' && [ "$#" -lt 5 ] || length=$5
cpu='8' && [ "$#" -lt 6 ] || cpu=$6
trimmomatic $mode -threads $cpu $in $out SLIDINGWINDOW:$window:$quality MINLEN:$length
