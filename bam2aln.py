#!/usr/bin/env python3

import sys
import re
import pysam
from argparse import ArgumentParser
from collections import defaultdict


cli = ArgumentParser(description='Convert BAM to Clustal style '
                     'alignments')
cli.add_argument('-b', '--bam', help='BAM file (default=STDIN)')
cli.add_argument('-r', '--region', help='reference region, e.g. chr1:7+81',
                 default='_:0+0')
cli.add_argument('-c', '--consensus', help='print consensus only',
                 action='store_true')
cli.add_argument('-t', '--threshold', help='base must have at least this '
                 'fraction of total coverage (default=0.6)', type=float,
                 default=0.6)
cli.add_argument('-m', '--min-reads', help='minimal read coverage '
                 '(default=3)', type=int, default=3)
args = cli.parse_args()


TR = {' ': ' ', 'N': 'N', '-': '-', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}


def revcomp(seq):
    return ''.join([TR[x.upper()] for x in seq[::-1]])


m = re.match('^(.*):(\d+)(-|\+)(\d+)$', args.region)
if not m:
    sys.stderr.write('ERROR: region "%s" is not valid!\n' % args.region)
    exit()
chro = m.group(1) if m.group(1) != '_' else None
start = int(m.group(2)) if m.group(2) != '0' else None
strand = m.group(3)
end = int(m.group(4)) if m.group(4) != '0' else None
theta = args.threshold

sam = pysam.AlignmentFile(args.bam, 'rb')

cons = ''
qual, cov = [], []
inserts = defaultdict(int)
pile = sam.pileup(chro, start, end)
for col in pile:
    if start is not None and (col.pos < start or end < col.pos):
        continue
    freq = defaultdict(int)
    # sys.stderr.write('cov @ %s: %s\n' % (col.pos, col.n))
    for read in col.pileups:
        base = ' '
        pos = read.query_position
        if read.is_del or read.is_refskip:
            base = '-'
        elif read.indel:
            base = read.alignment.query_sequence[pos:(pos+read.indel)]
        else:
            base = read.alignment.query_sequence[pos]
        freq[base] += 1

    if len(freq) < 1:
        cons += 'N'
        qual.extend([0])
        cov.extend([0])
        continue
    best = ['N', 0]
    N = sum(freq.values())
    if args.min_reads <= N:
        for c, n in freq.items():
            if best[1] < n and theta < n/N:
                best = [c, n]
                if len(c) > 1:
                    inserts[col.pos] = len(c)-1
    cons += best[0]
    qual.extend(([best[1]/N] if N else [0]) * len(best[0]))
    cov.extend([N] * len(best[0]))
    q = [33+int(40*(best[1]/N-0.25)*4/3)] if N > 120 else [33+int(N/3*(best[1]/N-0.25)*4/3)]
    if len(best[0]) > 1:
        q *= len(best[0])
    # qual.extend(q)

if strand == '-':
    cons = revcomp(cons)
    qual = qual[::-1]
    cov = cov[::-1]

print('@consensus_quality\n%s\n+\n%s\n@consensus_coverage\n%s\n+\n%s' %
      (cons, ''.join([chr(33+int(40*(q-0.25)*4/3)) for q in qual]),
       cons, ''.join([chr(33+int(c/3)) if c < 120 else 'H' for c in cov])))

# reads = []
# offset = 0
# columns = defaultdict(list)
# for a in sam.fetch(chro, start, end):
#     reads.append(a)
# for i in range(start, end):
#     for r in reads:
#         if start <= r.reference_start and r.reference_start+r.template_length <= end:
#             print(r.cigar)
sam.close()
