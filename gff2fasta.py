#!/usr/bin/python3

from Bio import SeqIO
import re
import sys

fasta_name = ''
gff_name = ''
padding = 0
maxlen = 999999
typefilter = None


def usage():
    u = """
gff2fasta - extraction of sequences from FASTA genomes with GFF annotations
USAGE: gff2fasta -f FASTA [-g GFF] [-p PADDING] [-l MAXLENGTH] [-t TYPES]
EXAMPLE: gff2fasta -f genome.fna -p 50 -l 1000 -t gene,ORF < annot.gff
"""
    sys.stderr.write(u)
    sys.stderr.flush()


if len(sys.argv) == 1:
    usage()
    quit()

i = 1
while len(sys.argv) > i:
    if sys.argv[i].startswith('-g'):
        gff_name = sys.argv[i+1]
    elif sys.argv[i].startswith('-f'):
        fasta_name = sys.argv[i+1]
    elif sys.argv[i].startswith('-t'):
        if ',' in sys.argv[i+1]:
            typefilter = set(sys.argv[i+1].split(','))
        else:
            typefilter = [sys.argv[i+1]]
    elif sys.argv[i].startswith('-p'):
        padding = int(sys.argv[i+1])
    elif sys.argv[i].startswith('-l'):
        maxlen = int(sys.argv[i+1])
    elif sys.argv[i].startswith('-h'):
        usage()
        quit()
    else:
        sys.stderr.write("unknown argument '%s'\n" % sys.argv[i])
        usage()
        quit()
    i += 2


fasta_seqs = SeqIO.to_dict(SeqIO.parse(fasta_name, 'fasta'))
p_id = re.compile(('^(\S+)\s+\S+\s+(\S+)\s+'        # chromosome, _, type
                   '(\d+)\s+(\d+)\s+\S+\s+([\+-])'  # from, to, _, strand
                   '\s+\S+\s+ID=([^;]+);(.+)$'))    # _, ID, rest
p = re.compile('^(\S+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+([\+-]).*$')

gff = sys.stdin
if gff_name != '':
    gff = open(gff_name, 'r')
else:
    sys.stderr.write("Waiting for GFF-formatted input:\n")

for line in gff:
    m = p.match(line)
    if m:
        chro = m.group(1)
        feature = m.group(2)
        start = int(m.group(3))
        end = int(m.group(4))
        strand = m.group(5)
        ident = '%s|%d-%d' % (chro, start, end)
        info = ''
        m = re.search(p_id, line)
        if m:
            ident = m.group(6)
            info = m.group(7)

        if end < start:
            start, end = (end, start)
        start -= padding
        end += padding
        if start < 1:
            start = 1
        if end > len(str(fasta_seqs[chro].seq)):
            end = len(str(fasta_seqs[chro].seq))

        if typefilter is None or feature in typefilter:
            seq = fasta_seqs[chro].seq[start:end]
            if strand is '-':
                seq = seq.reverse_complement()

            if abs(start-end) < maxlen:
                print(">%s|%d-%d|%s|%s\n%s\n" %
                      (fasta_seqs[chro].name,
                       start, end, strand,
                       ident, seq))

if gff_name != '':
    gff.close()
