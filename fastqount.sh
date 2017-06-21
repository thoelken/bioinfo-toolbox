#!/bin/bash
cmd="cat"
[[ $1 = *.gz ]] || [[ $1 = *.bz2 ]] && cmd="zcat"
echo "$(( $($cmd $1 | wc -l) / 4 ))"
