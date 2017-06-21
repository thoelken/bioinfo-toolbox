#!/usr/bin/env sh
makeblastdb -dbtype nucl -in $2 > /dev/null
blastn -db $2 -query $1 -outfmt 6 $3
