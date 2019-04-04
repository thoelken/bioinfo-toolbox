#!/usr/bin/env python3
from argparse import ArgumentParser as AP
from collections import defaultdict
import sys
import re


if sys.stdin.isatty():
    sys.stderr.write('waiting for SAM formatted input...\n')
    sys.stderr.flush()


cli = AP('Generates consensus sequence of SAM alignment')
cli.add_argument('-v', '--verbose', action='store_true', help='show debugging info')
cli.add_argument('-r', '--reads', type=int, default=3,
                 help='minimal number of reads per position [default: 3]')
cli.add_argument('-t', '--threshold', type=float, default=0.6,
                 help='minimal fraction of reads [default: 0.6]')
cli.add_argument('-F', '--formatting', choices=['fastq', 'fasta', 'clustal', 'seq'], default='seq',
                 help='output format: fastq, fasta, clustal, seq [=default]')
cli.add_argument('-R', '--region', help='process only defined region e.g. "chr1:1-200"')
args = cli.parse_args()

if args.region:
    chro, positions = args.region.split('|')
    start, end = [int(x) for x in positions.replace('+', '-').split('-')]
    antisense = '-' in positions

TR = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N', '-': '-'}


def reverse_complement(sequence):
    return ''.join([TR[x] for x in sequence[::-1]])


align = args.formatting.lower() in ['clustal', 'fasta']
if args.formatting.lower() == 'clustal':
    print('CLUSTAL W non-standard clustalW-like output\n')

data = {}
line_count = 0
no_header = True
for l in sys.stdin:
    line_count += 1
    if l.startswith('@'):
        if l.startswith('@SQ'):
            tags = {k: v for k, v in [tag.split(':') for tag in l.rstrip().split('\t')[1:]]}
            if args.region:
                if tags['SN'] == chro:
                    data[chro] = defaultdict(lambda: defaultdict(int))
                    data[chro]['start'] = start
                    data[chro]['end'] = end
            else:
                data[tags['SN']] = defaultdict(lambda: defaultdict(int))
                data[tags['SN']]['start'] = 1
                data[tags['SN']]['end'] = int(tags['LN'])
            no_header = False
        continue
    if no_header:
        raise Exception('SAM HEADER MISSING: use "samtools view -h" for input')
    c = l.rstrip().split('\t')
    if len(c) < 9:
        raise Exception('BAD SAM FORMAT: less than 9 columns in line %d' % line_count)
    if c[2] in data:
        s, e = data[c[2]]['start'], data[c[2]]['end']
        cigar = ''.join([x[1]*int(x[0]) for x in re.findall('(\d+)(\w)', c[5])])
        matches = ''.join([x for x in cigar if x not in ['S', 'I']])
        pos = int(c[3])
        if args.region and (pos+len(matches) < start or pos > end):
            continue
        if align:
            head = 'R%d %s\t' % (line_count, c[0])
            if args.formatting.lower() == 'fasta':
                head = '>%s\n' % c[0]
            sys.stdout.write('%s%s' % (head, '-'*(pos-s)))
        insert = ''
        for cig, seq in zip(cigar, c[9]):
            if cig == 'S':
                continue
            if cig == 'I':
                insert += seq
                continue
            if insert:
                if align and s <= pos < pos+len(insert) <= e:
                    sys.stdout.write(insert)
                data[c[2]][pos][insert] += 1
                pos += 1
                insert = ''
            if cig == 'D':
                if align and s <= pos <= e:
                    sys.stdout.write('-')
                data[c[2]][pos]['-'] += 1
            else:
                if align and s <= pos <= e:
                    sys.stdout.write(seq)
                data[c[2]][pos][seq] += 1
            pos += 1
        if align:
            sys.stdout.write('-'*(e-pos) + '\n')


def qual_trans(q, c):
    return chr(33+int(40*(c/120*(q-0.25)*4/3 if c < 120 else (q-0.25)*4/3)))


for ref in data:
    if args.formatting.lower() == 'fastq':
        print('@%s:%d-%d' % (ref, data[ref]['start'], data[ref]['end']))
    elif args.formatting.lower() == 'fasta':
        print('>consensus')
    elif args.formatting.lower() == 'clustal':
        sys.stdout.write('consensus\t')
    qual, seq = '', ''
    for i in range(data[ref]['start'], data[ref]['end']+1):
        best, q = 'N', chr(33)
        if i in data[ref]:
            best = max(data[ref][i], key=data[ref][i].get)
            total = sum(data[ref][i].values())
            if data[ref][i][best] < args.reads or data[ref][i][best]/total < args.threshold:
                best = 'N'
            else:
                q = qual_trans(data[ref][i][best]/total, data[ref][i][best])
        qual += q
        seq += best
        # sys.stdout.write(best)
    if args.region and antisense:
        seq = reverse_complement(seq)
    sys.stdout.write(seq)
    sys.stdout.write('\n')
    if args.formatting.lower() == 'fastq':
        print('+\n%s' % qual)
    elif args.formatting.lower() == 'clustal':
        print('!quality\t%s' % qual)
