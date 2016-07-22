#!/bin/bash
bam=$(echo $1 | sed 's/.(sam|bam)$//')
samtools view -buSh - | samtools sort - $bam; samtools index $bam.bam
