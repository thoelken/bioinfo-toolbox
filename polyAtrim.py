#!/usr/bin/env python3
from argparse import ArgumentParser as AP
import sys
import re


cli = AP(description='simple poly(A) trimmer script for FASTA and FASTQ')
cli.add_argument('-N', '--nucleotide', default='A',
                 help='nucleotide or sequence that should be cut (default: A)')
cli.add_argument('-r', '--repetitions', type=int, default=10,
                 help=('minimal number of repitions (default: 10)'
                       'N=0 trims only occurrences at the end)'))


args = cli.parse_args()
PATTERN = re.compile('%s{%d}.*$' % (args.nucleotide, args.repetitions))
if args.repetitions == 0:
    PATTERN = re.compile('%s+?$' % args.nucleotide)
count = 0
edited = 0


def read_fasta(l, stream):
    global count, edited
    print(l)
    l = stream.readline().rstrip()
    n = re.sub(PATTERN, '', l)
    print(n)
    count += 1
    if not len(l) == len(n):
        edited += 1


def read_fastq(l, stream):
    global count, edited
    print(l)
    l = stream.readline().rstrip()
    n = re.sub(PATTERN, '', l)
    print(n)
    print(stream.readline().rstrip())
    l = stream.readline().rstrip()
    count += 1
    if not len(l) == len(n):
        edited += 1
        l = l[:len(n)]
    print(l)


stream = sys.stdin
sys.stderr.write('waiting for FASTQ/FASTA input...')
sys.stderr.flush()
first = True
for l in stream:
    if first:
        sys.stderr.write('\r                                 \r')
        first = False
    l = l.rstrip()
    if l.startswith('>'):
        read_fasta(l, stream)
    elif l.startswith('@'):
        read_fastq(l, stream)
sys.stderr.write('  trimmed %d of %d sequences (%.2f%%)\n' % (edited, count, edited/count*100))
