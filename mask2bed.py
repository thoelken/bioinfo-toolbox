#!/usr/bin/env python3
from argparse import ArgumentParser as AP, FileType
import re

cli = AP()
cli.add_argument('-s', '--soft_masked', action='store_true',
                 help='Reference is softmasked (lower case "acgt")')
cli.add_argument('-f', '--fasta', type=FileType('r'), default='-', metavar='FILE',
                 help='Reference FASTA file [default: STDIN]')
args = cli.parse_args()

lowers = str.maketrans('acgt', 'nnnn')


def print_bed(head, seq, soft_masked=True):
    if soft_masked:
        seq = seq.translate(lowers)
    for m in re.finditer('n+', seq.lower()):
        print('%s\t%d\t%d' % (head, m.start(), m.end()))


head, seq = 'noRef', ''
for line in args.fasta:
    line = line.rstrip()
    if not line or line.startswith('!'):
        continue
    if line.startswith('>'):
        if seq:
            print_bed(head, seq, args.soft_masked)
        head = line.split()[0][1:]
        seq = ''
        continue
    seq += line
if seq:
    print_bed(head, seq, args.soft_masked)
