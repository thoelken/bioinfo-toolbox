sra=$1
name=$(echo $sra | sed 's/.sra//')
fasterq-dump $sra
if [ -e ${sra}_1.fastq ]; then
    cat ${sra}_1.fastq | gzip > $name.1.fq.gz
    cat ${sra}_2.fastq | gzip > $name.2.fq.gz
    rm ${sra}_1.fastq ${sra}_2.fastq
else
    cat $sra.fastq | gzip > $name.fq.gz
    rm $sra.fastq
fi
