#! /usr/bin/env bash
sra=$(echo $1 | sed 's/.sra//')
/home/clemens/bin/sratoolkit.2.5.2-ubuntu64/bin/sam-dump --header $1 | sam2bam $sra