#!/bin/bash
bam=$1; region=$2; ref=$3
gff=$(echo $2 | sed -r 's/^(.+):(.+)-(.+)$/\1\t.\tgene\t\2\t\3\t.\t+\t.\tID=ref;\n/')
#printf "$gff"
clustalw_many2one <(printf "$gff" | gff2fasta -f $ref) <(bam2fna $bam $region)
