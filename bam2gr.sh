#!/bin/bash
samtools view -b -F16 $1 $2 \
	| samtools depth /dev/stdin \
	| awk '{print $2"\t"$3}' \
	| cat - <(samtools view -b -f16 $1 $2 \
		| samtools depth /dev/stdin \
		| awk '{print $2"\t-"$3}')
