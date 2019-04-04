#!/bin/bash

param=""
if [[ $1 == \-* ]]; then
    param="$1"; shift;
fi

deep_join() {
    if [ $# -eq 1 ]; then
        join -t$'\t' $param - "$1";
    else
        f1=$1; shift;
        join -t$'\t' $param - "$f1" | deep_join "$@"
    fi
}

if [ $# -le 2 ]; then
    join -t$'\t' $param "$@";
else
    f1=$1; f2=$2; shift 2;
    join -t$'\t' $param "$f1" "$f2" | deep_join "$@";
fi
