#!/bin/bash

awk '{print $1"\tbed2gff\t"$4"\t"$2"\t"$3"\t"$5"\t"$6"\t.\tID=bed"NR";"}' $1
