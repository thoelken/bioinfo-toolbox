#!/bin/bash
samtools view -buSh - | samtools sort - > $1; samtools index $1
