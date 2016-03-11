#!/bin/bash
samtools view -b -F16 $1 \
	| samtools depth "${@:2}" /dev/stdin \
	| awk '{print $2"\t"$3}' \
	| cat - <(samtools view -b -f16 $1 \
		| samtools depth "${@:2}" /dev/stdin \
		| awk '{print $2"\t-"$3}')
