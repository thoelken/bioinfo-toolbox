#!/usr/bin/python3

from Bio import SeqIO
import re
import sys


p_id = re.compile(('^(\S+)\s+\S+\s+(\S+)\s+'        # chromosome, _, type
                   '(\d+)\s+(\d+)\s+\S+\s+([\+-])'  # from, to, _, strand
                   '\s+\S+\s+ID=([^;]+);(.+)$'))    # _, ID, rest
p = re.compile('^(\S+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+([\+-]).*$')
gip = re.compile('gi\|\d+\|ref\|(NC_\d+\.\d+)\|\s*.*')


def usage():
    u = """
gff2fasta - extraction of sequences from FASTA genomes with GFF annotations
USAGE: gff2fasta -f FASTA [-g GFF] [-p PADDING] [-l MAXLENGTH] [-t TYPES]
EXAMPLE: gff2fasta -f genome.fna -p 50 -l 1000 -t gene,ORF < annot.gff
"""
    sys.stderr.write(u)
    sys.stderr.flush()


def parse_gff(fasta_name, gff=sys.stdin, padding=0, maxlen=99999999, typefilter=None):
    fasta_seqs = SeqIO.to_dict(SeqIO.parse(fasta_name, 'fasta'))

    # sanitizing chromosome names
    tmp = dict()
    changed = 0
    for f in fasta_seqs:
        m = gip.match(f)
        if m:
            tmp[m.group(1)] = fasta_seqs[f]
            sys.stderr.write('found GI for refseq %s\n' % f)
            changed = 1
    if changed == 1:
        fasta_seqs = tmp

    for line in gff:
        m = p.match(line)
        if m:
            chro = m.group(1)
            feature = m.group(2)
            start = int(m.group(3))
            end = int(m.group(4))
            strand = m.group(5)
            ident = '%s|%d-%d' % (chro, start, end)
            m = re.search(p_id, line)
            if m:
                ident = m.group(6)
                # info = m.group(7)

            if end < start:
                start, end = (end, start)
            start -= padding
            end += padding
            if start < 1:
                start = 1
            if end > len(str(fasta_seqs[chro].seq)):
                end = len(str(fasta_seqs[chro].seq))

            if typefilter is None or feature in typefilter:
                seq = fasta_seqs[chro].seq[start-1:end]
                if strand is '-':
                    seq = seq.reverse_complement()

                if abs(start-end) < maxlen:
                    yield (chro, # fasta_seqs[chro].name,
                           start, end, strand,
                           ident, seq)


def print_fasta(seq_tuple):
    print(">%s|%d-%d|%s|%s\n%s\n" % seq_tuple)


def main():
    if len(sys.argv) == 1:
        usage()
        quit()

    fasta_name = ''
    gff_name = ''
    padding = 0
    maxlen = 99999999
    typefilter = None

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

    gff = sys.stdin
    if gff_name != '':
        with open(gff_name, 'r') as gff:
            for g in parse_gff(fasta_name, gff, padding, maxlen, typefilter):
                print_fasta((g))
    else:
        sys.stderr.write("Waiting for GFF-formatted input:\n")
        for g in parse_gff(fasta_name, sys.stdin, padding, maxlen, typefilter):
            print_fasta((g))


if __name__ == "__main__":
    main()
