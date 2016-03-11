#!/bin/bash
samtools view $1 $2 | awk -F $'\t' '{print ">"NR" "$1"\n"$10 }'
#'{ print "echo \">"$1"\""; rev=""; if(and($2, 0x10)){rev=" | tr ACGT TGCA | rev"} print "echo "$10""rev }' | bash
