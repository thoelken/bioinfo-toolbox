#!/bin/bash
# FASTQ subsampling based on http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/subsampling_reads.pdf
set -euo pipefail
if [[ $# -lt 2 || $1 == "-h" ]]; then
    echo "USAGE: $(basename $0) K FILE [FILE2]"
    printf "\tK:\tnumber of reads to sample (not precise for BAM and SAM files)\n"
    printf "\tFILE:\tfile with all reads (first file if paired end) [bam,sam,fq,fq.gz]\n"
    printf "\tFILE:\tfile for second read in pair (first file must be same type), [fq,fq.gz]\n"
    echo "EXAMPLE: $(basename $0) 1e6 alignment.bam"
    echo "EXAMPLE: $(basename $0) 1000000 reads.1.fq.gz reads.2.fq.gz"
    exit
fi
n=$(echo $2 | sed -r 's/\.(fastq|fq|sam|bam)(\.gz)?/.sub/')
if [[ $2 == *".sam" ]]; then
    frac=$(echo "$1 / $(samtools view -Sc $2)" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l )
    paired=$(samtools view -Scf 1 $2)
    if [[ $paired -eq 0 ]]; then
        samtools view -\@$(nproc) -Sbs $frac $2 | samtools fastq -\@$(nproc) - > $n.fq
    else
        samtools view -\@$(nproc) -Sbs $frac $2 | samtools fastq -\@$(nproc) -n -1 $n.1.fq.gz -2 $n.2.fq.gz -
    fi
elif [[ $2 == *".bam" ]]; then
    frac=$(echo "$1 / $(samtools view -c $2)" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l )
    paired=$(samtools view -c -f 1 $2)
    if [[ $paired -eq 0 ]]; then
        samtools view -\@$(nproc) -bs $frac $2 | samtools fastq -\@$(nproc) - | gzip > $n.fq.gz
    else
        frac=$(echo "$frac * 2" | bc -l)
        samtools view -\@$(nproc) -bs $frac $2 | samtools fastq -\@$(nproc) -n -1 $n.1.fq.gz -2 $n.2.fq.gz -
    fi
elif [[ $2 == *".gz" ]]; then
    if [[ $# -eq 3 ]]; then
        m=$(echo $3 | sed -r 's/\.(fastq|fq|sam|bam)(\.gz)?/.sub/')
        paste <(zcat $2) <(zcat $3) |\
            awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' |\
            awk -v k=$1 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |\
            awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 | "gzip > '$n'.fq.gz";print $2"\n"$4"\n"$6"\n"$8 | "gzip > '$m'.fq.gz"}'
    else
        zcat $2 | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' |\
            awk -v k=$1 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |\
            awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4}' | gzip > $n.fq.gz
    fi
elif [[ $2 == *".fq" || $1 == *".fastq" ]]; then
    if [[ $# -eq 3 ]]; then
        m=$(echo $3 | sed -r 's/\.(fastq|fq|sam|bam)(\.gz)?/.sub/')
        paste $2 $3 |\
            awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' |\
            awk -v k=$1 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |\
            awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "'$n'.fq";print $2"\n"$4"\n"$6"\n"$8 > "'$m'.fq"}'
    else
        cat $2 | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' |\
            awk -v k=$1 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |\
            awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4}' > $n.fq
    fi
fi
