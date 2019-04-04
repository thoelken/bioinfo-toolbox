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
cli.add_argument('-a', '--align', help='print line by line pseudo alignment', action='store_true')
cli.add_argument('-c', '--consensus', help='print consensus only', action='store_true')
cli.add_argument('-t', '--threshold', help='base must have at least this '
                 'fraction of total coverage (default=0.6)', type=float,
                 default=0.6)
cli.add_argument('-m', '--min-reads', help='minimal read coverage '
                 '(default=3)', type=int, default=3)
args = cli.parse_args()


TR = {' ': ' ', 'N': 'N', '-': '-', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
CIGAR = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']


def revcomp(seq):
    return ''.join([TR[x.upper()] for x in seq[::-1]])


def make_aln(aln, start=0, end=0, antisense=False):
    line = ''
    pos = aln.reference_start
    seq = aln.query_sequence
    pairs = {r: q for q, r in aln.get_aligned_pairs()}
    cigar = ''.join([CIGAR[a]*n for a, n in aln.cigartuples])
    cigar = {pos+i: a for i, a in enumerate(cigar)}
    for i in range(start, end+1):
        if i in pairs and pairs[i] is not None and i in cigar and cigar[i] != 'S':
            line += seq[pairs[i]-1]
        else:
            line += '-'
    if antisense:
        line = revcomp(line)
    return line
    # if pos+len(seq) <= start or end <= pos:
    #     return None
    # if start < pos:
    #     line += '-'*(pos-start)
    # if end < pos+len(seq):
    #     seq = seq[:len(seq)-(pos+len(seq)-end)]
    # if pos < start:
    #     seq = seq[start-pos:]
    # line += seq
    # if len(line) < end-start:
    #     line += '-'*(end-start-len(line))
    # if antisense:
    #     line = revcomp(line)
    # return line


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

if args.align:
    print('CLUSTAL W non-standard clustalW-like output\n')
    reads = {}
    for col in sam.pileup(chro, start, end):
        for read in col.pileups:
            if read.alignment.query_name not in reads:
                reads[read.alignment.query_name] = 1
                line = make_aln(read.alignment, start, end, strand == '-')
                if line is not None:
                    print(read.alignment.query_name + '\t' + line)

if args.consensus:
    cons = ''
    qual, cov = [], []
    inserts = defaultdict(int)
    last = start
    for col in sam.pileup(chro, start, end):
        if start is not None and (col.pos < start or end < col.pos):
            continue  # skip columns out of range
        if last < col.pos-1:  # skipped some columns without reads
            cons += 'N' * (col.pos-last)
            qual.extend([0.25] * (col.pos-last))
            cov.extend([0] * (col.pos-last))
        last = col.pos
        freq = defaultdict(int)
        for read in col.pileups:
            base = ' '
            pos = read.query_position
            if read.is_del or read.is_refskip:  # deletion in query
                base = '-'
            elif read.indel:  # insertion in query
                base = read.alignment.query_sequence[pos:(pos+read.indel)]
            else:
                base = read.alignment.query_sequence[pos]
            freq[base] += 1

        N = sum(freq.values()) if freq else 0
        if N < args.min_reads:
            cons += 'N'
            qual.extend([0.25])
            cov.extend([N])
            continue
        best = ['N', 1]
        for c, n in freq.items():
            if best[1] < n:
                best = [c, n]
                if len(c) > 1:
                    inserts[col.pos] = len(c)-1
        cons += best[0] if theta < best[1]/N else 'N'
        qual.extend([best[1]/N] * len(best[0]))
        cov.extend([N] * len(best[0]))
    if last < end:
        cons += 'N' * (end-last)
        qual.extend([0.25] * (end-last))
        cov.extend([0] * (end-last))

    if strand == '-':
        cons = revcomp(cons)
        qual = qual[::-1]
        cov = cov[::-1]

    def qual_trans(q, c):
        return chr(33+int(40*(c/120*(q-0.25)*4/3 if c < 120 else (q-0.25)*4/3)))

    quality = ''.join([qual_trans(q, c) for q, c in zip(qual, cov)])
    if args.align:
        print('consensus\t%s\nidentity\t%s\ncoverage\t%s' % (cons, quality, quality))
        #      (cons, ''.join([chr(33+int(40*(q-0.25)*4/3)) for q in qual]),
        #        ''.join([chr(33+int(c/3)) if c < 120 else 'I' for c in cov])))
    else:
        print('@consensus %s\n%s\n+\n%s' % (args.region, cons, quality))

sam.close()
