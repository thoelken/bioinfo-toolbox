#!/usr/bin/env python3
import sys
from argparse import ArgumentParser as AP
import logging as log

cli = AP(description='Change annotation positions based on multiple sequence alignment')
cli.add_argument('-v', '--verbose', action='store_true', help='output additional infos to STDERR')
cli.add_argument('-t', '--tags', nargs='+', help='output only certain tags from last GFF column ' +
                 '[default: all]')
cli.add_argument('-o', '--outprefix', help='print out all features into files with this prefix')
cli.add_argument('-m', '--merge', action='store_true', help='merge overlapping features')
args = cli.parse_args()
log.basicConfig(level=log.DEBUG if args.verbose else log.WARN, format='[%(levelname)s] %(message)s')


def parse_gff(stream):
    """add genes from GFF into global variable by position and
    return genome specific subset"""
    genes = []
    for l in stream:
        l = l.rstrip()
        if not l or l.startswith('#'):
            continue
        c = l.split('\t')
        if len(c) < 9:
            continue
        if not args.tags or c[2] in args.tags:
            genes.append(c)
    log.info('parsed %d unique features' % len(genes))
    return genes


genes = parse_gff(sys.stdin)
if args.offset is None:
    offset = min([s for chro in genes for s in genes[chro]])
