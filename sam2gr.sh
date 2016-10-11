#!/bin/bash
samtools view -Sb -F16 $1 $2 \
	| samtools depth /dev/stdin \
	| awk '{print $2"\t"$3}' \
	| cat - <(samtools view -Sb -f16 $1 $2 \
		| samtools depth /dev/stdin \
		| awk '{print $2"\t-"$3}')
