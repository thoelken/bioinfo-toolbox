#!/usr/bin/env python3
import sys

n, t, o = 0, 0, 0
gene, trans = '', ''
for l in sys.stdin:
    n += 1
    if not l or l.startswith('#'):
        continue
    c = l.rstrip().split("\t")
    if len(c) != 9:
        sys.stderr.write('WARNING: not 9 columns in line: %d' % n)
        continue
    atr = {x.lstrip().split(' ', 1)[0]: x.lstrip().split(' ', 1)[1][1:-1] for x in c[8][:-1].split(';')}
    ident = atr.pop('gene_id')
    if c[2] == 'gene':
        gene = ident
        trans = ''
        t = 0
        o = 0
    elif c[2] in ['mRNA', 'transcript']:
        atr['parent'] = gene
        t += 1
        ident += '_t%d' % t
        trans = ident
    else:
        atr['parent'] = trans if trans else gene
        o += 1
        ident += '_o%d' % o
    atr['ID'] = ident
    c[8] = ';'.join(['%s=%s' % (k, v) for k, v in atr.items()])
    print('\t'.join(c))
