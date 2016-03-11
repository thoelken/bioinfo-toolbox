#!/usr/bin/env python3

import sys
import re
import pysam
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description='Convert BAM to Clustal style '
                                 'alignments')
parser.add_argument('-b', '--bam', help='BAM file (default=STDIN)')
parser.add_argument('-r', '--region', help='reference region, e.g. chr1:7+81',
                    default='_:0+0')
parser.add_argument('-c', '--consensus', help='print consensus only',
                    action='store_true')
parser.add_argument('-t', '--threshold', help='base must have at least this '
                    'fraction of total coverage (default=0.6)', type=float,
                    default=0.6)
parser.add_argument('-m', '--min-reads', help='minimal read coverage '
                    '(default=3)', type=int, default=3)


args = parser.parse_args()

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
pile = sam.pileup(chro, start, end)
for col in pile:
    if start is not None and (col.pos < start or end < col.pos):
        continue
    freq = defaultdict(int)
    inserts = defaultdict(int)
    # sys.stderr.write('cov @ %s: %s\n' % (col.pos, col.n))
    for read in col.pileups:
        if read.is_del or read.is_refskip:
            freq['-'] += 1
        elif read.indel:
            pos = read.query_position
            freq[read.alignment.query_sequence[pos:(pos+read.indel)]] += 1
        else:
            freq[read.alignment.query_sequence[read.query_position]] += 1
    if theta < sum(inserts.values()) / len(col.pileups):
        m = max(inserts.values())
        for i in inserts:
            if inserts[i] == m:
                sys.stdout.write(i)
                break

    best = ['X', 0]
    N = sum(freq.values())
    if args.min_reads <= N:
        for c, n in freq.items():
            if best[1] < n and theta < n/N:
                best = [c, n]

    sys.stdout.write(best[0])
sys.stdout.write('\n')
sam.close()
