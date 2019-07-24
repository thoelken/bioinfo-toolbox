#!/usr/bin/env python3
from argparse import ArgumentParser as AP
import sys
import re


cli = AP(description='simple poly(A) trimmer script for FASTA and FASTQ')
cli.add_argument('-N', '--nucleotide', default='A',
                 help='nucleotide or sequence that should be cut (default: A)')
cli.add_argument('-r', '--repetitions', type=int, default=8,
                 help=('minimal number of repetitions '
                       '(N=0 trims only occurrences at the end; default: 8)'))
cli.add_argument('-l', '--min_length', type=int, default=20, metavar='N',
                 help='discard reads with length < N after trimming (default: 20)')
# cli.add_argument('files', nargs='*', help='input FASTQ file(s) (default: read from STDIN)')


args = cli.parse_args()
PATTERN = re.compile('%s{%d}.*$' % (args.nucleotide, args.repetitions))
if args.repetitions == 0:
    PATTERN = re.compile('%s+?$' % args.nucleotide)
count, edited, removed = 0, 0, 0


def read_fasta(l, stream):
    global count, edited, removed
    h = l
    l = stream.readline().rstrip()
    s = re.sub(PATTERN, '', l)
    if len(s) < args.min_length:
        removed += 1
        return
    print('%s\n%s' % (h, s))
    count += 1
    if not len(l) == len(s):
        edited += 1


def read_fastq(l, stream):
    global count, edited, removed
    h = l
    l = stream.readline().rstrip()
    s = re.sub(PATTERN, '', l)
    c = stream.readline().rstrip()
    q = stream.readline().rstrip()
    count += 1
    if len(s) < args.min_length:
        removed += 1
        return
    if not len(q) == len(s):
        edited += 1
        q = q[:len(s)]
    print('%s\n%s\n%s\n%s' % (h, s, c, q))


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
sys.stderr.write('total reads: %d\nremoved: %d (%.2f%%)\ntrimmed %d (%.2f%%)\n' %
                 (count, removed, removed/count*100, edited, edited/count*100))
