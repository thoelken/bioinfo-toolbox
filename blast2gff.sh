#!/usr/bin/env sh
# sed -rn 's/^\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+([0-9]+)\s+([0-9]+)\s+(\S+)\s+(\S+).*$/\1\tBLAST\thit\t\2\t\3\t.\t+\t.\tID=id0;Evalue=\4;bitscore=\5;/p'
awk '{strand="+"; if($9>$10){strand="-"; tmp=$9; $9=$10; $10=tmp}; print $2 "\tBLAST\thit\t" $9 "\t" $10 "\t.\t" strand "\t.\tID=id" NR ";Evalue=" $11 ";Bitscore=" $12 ";Identity=" $3 ";"}'
