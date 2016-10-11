#!/usr/bin/env sh
sed -rn 's/^.*\(([0-9]+)\)\s+!\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+([0-9]+)\s+([0-9]+)\s+(-|\+)\s+.*$/\4\tinfernal\tprediction\t\5\t\6\t.\t\7\t.\tID=pred\1;Evalue=\2;score=\3;/p' | awk 'BEGIN{OFS="\t"} {if($4>$5) {tmp=$4;$4=$5;$5=tmp}; print $0}'
