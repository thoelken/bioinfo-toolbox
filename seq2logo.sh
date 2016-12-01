#!/bin/bash
if [ -t 0 ]; then
    echo "USAGE: seq2logo 20 1000 <reads.fastq >logo.ps    # for length of 20 in first 1000 sequences"
    exit
fi
[ -n "$1" ] && l=$1 || l=20
[ -n "$2" ] && n=$(($2 * 2)) || n=2000
awk '
/^[>@]/ {
    fastq=sub(/@/, ">");
    print $0;
    getline;
    if(length($0) > '$l') {
        $0 = substr($0, 0, '$l');
        print $0
    } else {
        printf $0;
        for (i=1; i+length($0) <= '$l'; i++)
            printf "-";
        printf "\n"
    } if(fastq==1) {
        getline;
        getline
    }
}' | head -n$n | weblogo
