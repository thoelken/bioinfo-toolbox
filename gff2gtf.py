#!/usr/bin/env python3

import sys
from argparse import ArgumentParser as AP
from collections import defaultdict, deque

cli = AP(description='converts GFF3 to GTF files with gene-transcript-exon relations')
cli.add_argument('-c', '--clean', action='store_true', help='clean up duplicates')
args = cli.parse_args()

lines = 0
loc2gene = {}
id2line = {}
relations = defaultdict(list)
last_gene, last_transcript = '', ''
for l in sys.stdin:
    lines += 1
    l = l.rstrip()
    if l.startswith('#'):
        print(l)
        continue

    # parsing columns and attributes
    c = l.split('\t')
    if len(c) < 9:
        sys.stderr.write('[ERROR] Wrong format in GFF file ( <9 columns) @ line %d\n' % lines)
        sys.exit(1)
    chro, typ, start, end, strand, attr = c[0], c[2], c[3], c[4], c[6], c[8]
    loc = chro+'|'+typ+'|'+start+strand+end
    attr = [a for a in attr.split(';')]
    tmp = {}
    for a in attr:
        k, v = a.split('=')
        tmp[k] = v
    attr = tmp
    ident = attr['ID']

    # processing gene model
    if 'ID' not in attr:
        sys.stderr.write('[ERROR] Feature without ID @ line %d\n' % lines)
        sys.exit(1)

    attr['gene_id'] = last_gene
    attr['transcript_id'] = last_transcript
    id2line[ident] = '\t'.join(c[0:7])+'\t'+' '.join([k+v for k, v in attr])
    if typ == 'gene':
        loc2gene[loc] = ident
        last_gene = ident
    elif typ in ['transcript', 'mRNA']:
        last_transcript = ident

    if 'Parent' in attr:
        relations[attr['Parent']].append(ident)
    elif 'Gene' in attr:
        relations[attr['Gene']].append(ident)


for l, g in sorted(loc2gene):
    print(id2line(g))
    r = deque(relations[g])
    while r:
        c = r.popleft()
        print(id2line[c])
        r.extendleft(relations[c])
