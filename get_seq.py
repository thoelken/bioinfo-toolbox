#!/usr/bin/env python3
import re
import sys
from argparse import FileType, ArgumentParser as AP

cli = AP('extracts sequences from FASTA based on GFF annotations')
cli.add_argument('-f', '--fasta', help='FASTA reference file', required=True)
cli.add_argument('-a', '--annotation', help='annotation file (GFF,GTF,BED)',
                 type=FileType('r'), default=sys.stdin)
cli.add_argument('-p', '--padding', help='extracted sequences are padded'
                 'by N nucleotides', type=int, default=0, metavar='N')
cli.add_argument('--first', help='first N bases of each sequence',
                 type=int, default=0, metavar='N')
cli.add_argument('-t', '--typeof', help='filter by type of feature in in the 3rd column of GFF',
                 nargs='*', default=[])
cli.add_argument('-s', '--min_score', help='filter by score', default=0, metavar='F', type=float)
cli.add_argument('-i', '--ident', help='identity field in last column',
                 default='ID')
cli.add_argument('-r', '--region', metavar='CHR:START-END',
                 help='extract sequence from START to END or region '
                 '(higher number first for anti-sense strand)')


ID_PATTERN = re.compile('^(\S+)\s+\S+\s+(\S+)\s+'         # chro, _, type
                        '(\d+)\s+(\d+)\s+\S+\s+([\+-])'   # from, to, _, strand
                        '\s+\S+\s+ID=([^;]+);?(.+)$')     # _, ID, rest
GFF_PATTERN = re.compile('^(\S+)\t+[^\t]+\s+(\S+)\s+(\d+)\s+'
                         '(\d+)\s+\S+\s+([\+-]).*$')
GENOME_PATTERN = re.compile('gi\|\d+\|ref\|(NC_\d+\.\d+)\|\s*.*')


def parse_gff_line(line):
    if line.startswith('#'):
        return None
    m = GFF_PATTERN.match(line)
    if m:
        chro = m.group(1)
        feature = m.group(2)
        start = int(m.group(3))
        end = int(m.group(4))
        strand = m.group(5)
        ident = '%s|%d-%d' % (chro, start, end)
        m = re.search(ID_PATTERN, line)
        if m:
            ident = m.group(6)
        if end < start:
            start, end = (end, start)
        return (chro, feature, start-1, end, strand, ident)


def parse_fasta(fasta):
    ret = {}
    head, seq = '', ''
    with open(fasta) as lines:
        for li in lines:
            if not li or li[0] in '!#':
                continue
            li = li.rstrip()
            if li.startswith('>'):
                if seq:
                    ret[head] = seq
                head = li[1:].split(' ')[0]
                seq = ''
            else:
                seq += li
    if seq:
        ret[head] = seq
    return ret


TR = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def revcomp(s):
    return s.translate(TR)[::-1]


def parse_gff(refs, gff=sys.stdin, padding=0, maxlen=0, typefilter=[]):
    if gff == sys.stdin:
        sys.stderr.write("reading annotation input from STDIN...")
    for line in gff:
        sys.stderr.write('\r                                       \r')
        gff_line = parse_gff_line(line)
        if gff_line is None:
            continue
        (chro, feature, start, end, strand, ident) = gff_line
        start -= padding
        end += padding
        if start < 1:
            start = 1
        if end > len(refs[chro]):
            end = len(refs[chro])
        if not typefilter or feature in typefilter:
            seq = refs[chro][start:end]
            if strand is '-':
                seq = revcomp(seq)
            print(">%s|%d-%d|%s|%s\n%s" % (chro, start, end, strand, ident, seq))


if __name__ == "__main__":
    args = cli.parse_args()
    ID_PATTERN = re.compile('^(\S+)\s+\S+\s+(\S+)\s+'         # chro, _, type
                            '(\d+)\s+(\d+)\s+\S+\s+([\+-])'   # from, to, _, strand
                            '\s+\S+\s+.*%s=([^;]+);?(.+)$' %  # _, ID, rest
                            args.ident)

    refs = parse_fasta(args.fasta)
    parse_gff(refs, args.annotation, args.padding, args.typeof)
