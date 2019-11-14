#!/usr/bin/env python3
from argparse import ArgumentParser as AP, FileType
import sys

cli = AP()
cli.add_argument('-s', '--soft_masked', action='store_true',
                 help='Reference is softmasked (lower case "acgt")')
cli.add_argument('fasta', type=FileType('r'), default=sys.stdin, help='Reference FASTA file')
args = cli.parse_args()

lowers = str.maketrans('acgt', 'nnnn')

head, index, rep = 'noRef', 0, -1
for line in args.fasta:
    line = line.rstrip()
    if not line or line.startswith('!'):
        continue
    if line.startswith('>'):
        if 0 <= rep:
            print('%s\t%d\t%d' % (head, rep, index))
        head = line.split()[0][1:]
        index, rep = 0, -1
        continue
    if args.soft_masked:
        line = line.translate(lowers)
    for c in line.lower():
        if 0 <= rep and c != 'n':
            print('%s\t%d\t%d' % (head, rep, index))
            rep = -1
        elif rep < 0 and c == 'n':
            rep = index
        index += 1
if 0 <= rep:
    print('%s\t%d\t%d' % (head, rep, index))
