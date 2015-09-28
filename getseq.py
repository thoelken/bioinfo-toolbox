#!/usr/bin/python

from Bio import SeqIO
import sys

fasta = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

fasta_file = open(fasta)
fasta_seq = SeqIO.parse(fasta_file, 'fasta').next()
seq = fasta_seq.seq[start:end]
print(">%s %d-%d %s %s\n%s\n" %
      (fasta_seq.name, start, end, strand, ident, seq))
    # name, sequence = fasta.id, fasta.seq.tostring()

fasta_file.close()
