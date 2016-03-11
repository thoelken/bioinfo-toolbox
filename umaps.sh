#!/bin/bash

samtools view -F 0x904 $1 \
	| fgrep -w NH:i:1 \
	| cut -f1 \
	| sort \
	| uniq \
	| wc -l
