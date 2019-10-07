#!/usr/bin/env python3
import pysam
from argparse import ArgumentParser as AP
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

cli = AP(description='Output bins of the reference in SAM/BAM alignment')
cli.add_argument('-b', '--binsize', type=int, metavar='N', default=1000000,
                 help='window size for binning [default: 1000000]')
cli.add_argument('-f', '--format', choices=['gtf', 'gff', 'bed'], default='gtf',
                 help='output format in {gtf,gff,bed} [default: gtf]')
cli.add_argument('-o', '--overlap', type=float, metavar='F', default=0.0,
                 help='overlap between windows [default: 0.0]')
cli.add_argument('--skip_long_name', action='store_true',
                 help='skip chromosomes with longer names than 6 characters [default: false]')
cli.add_argument('--plus_strand_only', action='store_true',
                 help='print strand as + instead of . in GFF/GTF')
cli.add_argument('bam', default='-', help='SAM/BAM alignment as reference (use "-" for STDIN)')
args = cli.parse_args()

if 1 < args.overlap:  # in case the user specifies overlap in number of bases
    args.overlap = args.overlap / args.binsize


def print_bin(chro, length, binsize):
    subbin = (1-args.overlap)*binsize
    bins = [i for i in range(int(length/subbin))]
    for b in bins:
        if args.format == 'bed':
            print('%s\t%d\t%d' % (chro, b*subbin, b*subbin+binsize-1))
            continue
        idtag = '\tgene_id "%s-%d"; transcript_id "%s-%d";' % (chro, b, chro, b)
        if args.format == 'gff':
            idtag = '\tID=%s-%d;' % (chro, b)
        print('%s\tbin\texon\t%d\t%d\t0\t.\t.\t%s' % (chro, b*subbin+1, b*subbin+binsize, idtag))
    if not bins:  # sequence is shorter than binsize
        bins = [0]  # creates a bin from 0 to length
    if bins[-1]*subbin+binsize < length:  # add last partial bin till length
        if args.format == 'bed':
            print('%s\t%d\t%d' % (chro, (bins[-1]+1)*subbin, length-1))
            return
        idtag = 'gene_id "%s-%d"; transcript_id "%s-%d";' % (chro, len(bins), chro, len(bins))
        if args.format == 'gff':
            idtag = 'ID=%s-%d;' % (chro, len(bins))
        print('%s\tbin\texon\t%d\t%d\t0\t.\t.\t%s' % (chro, (bins[-1]+1)*subbin+1, length, idtag))


bam = pysam.AlignmentFile(args.bam)
for r in bam.references:
    print_bin(r, bam.get_reference_length(r), args.binsize)

# import sys
# for l in sys.stdin:
#     if not l.startswith('@'):  # seems to be end of header
#         break
#     if not l.startswith('@SQ'):  # parse only sequence length information
#         continue
#     c, bp = l.rstrip().split('\t')[1:3]
#     c = c[3:]
#     if args.skip_long_name and len(c) > 6:  # skip stuff like chr1_unplaced3298748
#         continue
#     bp = int(bp[3:])
#     print_bin(c, bp, args.binsize)
