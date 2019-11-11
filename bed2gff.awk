#!/usr/bin/awk -f
BEGIN{FS="\t"; OFS="\t"}
!/^#/{
    name = "bed"NR; score = 0; strand = "."
    if(NF>3) {name = $4}
    if(NF>4) {score = $5}
    if(NF>5) {strand = $6}
    print $1"\tbed2gff\tgene\t"($2+1)"\t"($3+1)"\t"score"\t"strand"\t.\tID="name";"
}
