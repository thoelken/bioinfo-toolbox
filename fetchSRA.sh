#!/bin/bash
s=$(echo $1 | sed 's/\..*$//')
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${s:0:3}/${s:0:6}/$s/$s.sra
