#!/usr/bin/env python3
import sys
import re
from argparse import ArgumentParser as AP


cli = AP(description='Offset SAM formated input to certain coordinates')
cli.add_argument('-o', '--offset', type=int, metavar='POS', default=0,
                 help='position to use as offset')
cli.add_argument('-r', '--reverse', action='store_true',
                 help='reverse mapping position (after offset)')
args = cli.parse_args()

for l in sys.stdin:
    if l.startswith('@'):
        print(l.rstrip())
    c = l.rstrip().split('\t')
    if len(c) < 11:
        continue
    pos = int(c[3])-args.offset
    if args.reverse:
        l = sum([int(x) for x in re.findall('(\d+)\w', c[5])])
        pos = -pos-l
        c[1] = str(int(c[1]) ^ 0x10)
        c[8] = str(-int(c[8]))
    c[3] = str(pos)
    print('\t'.join(c))
