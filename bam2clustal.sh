#!/bin/bash
samtools view $1 $2 | awk -F $'\t' '{print ">"$1"\n"$10}' | clustalo --pileup -i -
