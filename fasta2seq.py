#!/usr/bin/python

from Bio import SeqIO
import sys

fasta = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

fasta_file = open(fasta, 'rU')
fasta_seq = SeqIO.read(fasta_file, 'fasta')

seq = fasta_seq.seq[(start-1):end]
if start > end:
    seq = fasta_seq.seq[(end-1):start].reverse_complement()

print(">%s|%d-%d\n%s" %
      (fasta_seq.name, start, end, seq))
    # name, sequence = fasta.id, fasta.seq.tostring()

fasta_file.close()
