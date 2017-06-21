#!/usr/bin/env python3

import pysam
from argparse import ArgumentParser as AP


cli = AP(description='Creates SVG alignments with coverage from SAM/BAM input')
cli.add_argument('alignment', help='SAM or BAM file ("-" for SAM formated STDIN)')
cli.add_argument('-r', '--region', help='specified region e.g. "chr1:324-5343"')
args = cli.parse_args()

width = 5
height = 12

algn = pysam.AlignmentFile(args.alignment)

offset = 0

if args.region:
    ref, pos = args.region.split(':')
    start, end = [int(x) for x in pos.split('-')]
    algn = algn.fetch(ref, start, end)
    offset = start


print('<svg>')
row = 1
for l in algn:
    q = 0
    col = l.reference_start-offset
    for c, n in l.cigartuples:
        if c == 0:  # match
            for i in range(q, q+n):
                print('<text x="%d" y="%d">%s</text>"' % (col*width, row*height, l.query_sequence[i]))
                col += 1
            q += n
        elif c == 1:    # insert
            print('<text x="%d" y="%d">%s</text>"' % ((col+0.5)*width, row*height, 'v'))
            col += 1
            q += n
        elif c == 2:    # delete
            for _ in range(n):
                print('<text x="%d" y="%d">-</text>"' % (col*width, row*height))
                col += 1
    row += 1
print('</svg>')
