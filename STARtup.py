#!/usr/bin/env python3
# Fork from 'Anand M' star_aligner.py boiled down to minimal parameters

import os
import sys
from argparse import ArgumentParser as AP


cli = AP(description="Wrapper for the STAR mapping tool with quick and sane parameters")
cli.add_argument('reference', metavar='REF', help='path to genome FASTA')
cli.add_argument('reads', help='path/to/read1 [path/to/read2]', nargs='+', metavar='READS')

cli.add_argument('-t', '--threads', default=os.cpu_count(), type=int, metavar='N',
                 help='use N threads for mapping [%d]' % os.cpu_count())
cli.add_argument('-o', '--outprefix', type=str,
                 help='prefix of output files (including full or relative path) '
                 '[same as reads plus ".STAR."]')
cli.add_argument('-a', '--annotation', type=str,
                 help='GTF file with annotations')
cli.add_argument('-i', '--index', help='(re-)create STAR index in the specified folder')
cli.add_argument('--overwrite_index', action='store_true', help='overwrite existing STAR index')
cli.add_argument('-q', '--quant', action='store_true',
                 help='quantification based on annotation')
cli.add_argument('--sam', choices=['sam', 'bam', 'none'], default='bam', type=str,
                 metavar='sam|bam|none', help='output format of Alignment [bam]')
cli.add_argument('--unsort', action='store_true',
                 help='alignment is NOT sorted by coordinate (less RAM)')
cli.add_argument('-u', '--unmapped', action='store_true', help='write unmapped reads to extra file')
cli.add_argument('-m', '--multimap', default=10, type=int, metavar='N',
                 help='generate only alignments for reads that map at most N times')
cli.add_argument('-p', '--passthrough', default='',
                 help='parameters passed through to STAR (use quotes)')
args = cli.parse_args()

if len(args.reads) > 2:
    sys.stderr.write('ERROR: more than 2 read files given!\n')
    quit()

out = ''
if args.outprefix:
    out = args.outprefix
    if not out.endswith('.'):
        out += '.'
else:
    out = args.reads[0]
    if out.endswith('.gz'):
        out = out[0:-2]
    out = out[0:out.rfind('.', out.rfind('/'))] + '.STAR.'
outprefix = '--outFileNamePrefix ' + out + ' '

gzip = ''
r = args.reads[0]
if r.endswith('.gz') or r.endswith('zip'):
    gzip = '--readFilesCommand zcat '
elif r.endswith('.bz2'):
    gzip = '--readFilesCommand bzcat '

fmt = '--outSAMtype '
if args.sam.lower() == 'sam':
    fmt += 'SAM '
    args.unsort = True
elif args.sam.lower() == 'bam':
    fmt += 'BAM '
elif args.sam.lower() == 'none':
    fmt += 'None '
fmt += 'Unsorted ' if args.unsort else 'SortedByCoordinate '

gdir = args.index if args.index else args.reference + '.STAR'
if gdir.endswith('/'):
    gdir = gdir[:-1]

params = '--runThreadN %d --genomeDir %s ' % (args.threads, gdir)
params += args.passthrough + ' '
params += ('--sjdbGTFfile %s ' % args.annotation) if args.annotation else ''

if not os.path.isdir(gdir):
    os.mkdir(gdir)
if not os.path.exists(gdir + '/SAindex'):
    sys.stderr.write('\n[STARtup] Could not find index for "%s" ...\n' % args.reference)
    args.overwrite_index = True
if args.overwrite_index:
    cmd = ('STAR --runMode genomeGenerate --genomeFastaFiles %s %s' % (args.reference, params))
    sys.stderr.write('[STARtup] Creating new index "%s" with command:\n\n' % gdir)
    sys.stderr.write(cmd + '\n\n')
    if not os.system(cmd) == 0:
    	sys.stderr.write('ERROR: STAR index could not be generated!\n')
    	quit()

reads = ' '.join(args.reads)
params += gzip
params += outprefix
params += '--quantMode GeneCounts ' if args.quant else ''
params += ('--outFilterMultimapNmax %d ' % args.multimap) if args.multimap else ''
params += '--outReadsUnmapped Fastx ' if args.unmapped else ''
params += fmt
params += '--readFilesIn %s' % ' '.join(args.reads)
cmd = 'STAR %s' % params
sys.stderr.write('\n[STARtup] Running STAR with the following parameters:\n')
sys.stderr.write(cmd + '\n\n')
if not os.system(cmd) == 0:
    sys.stderr.write('ERROR: STAR mapping not successfull!\n')
    quit()

if args.sam.lower() == 'sam':
    if args.unsort:
        os.rename(out+'Aligned.out.sam', out+'sam')
    else:
        os.rename(out+'Aligned.sortedByCoord.out.sam', out+'sam')
if args.sam.lower() == 'bam' and not args.unsort:
    os.rename(out+'Aligned.sortedByCoord.out.bam', out+'bam')
    sys.stderr.write('\n[STARtup] Let me index that BAM for you ...\n')
    try:
        os.system('samtools index %sbam' % out)
    except OSError:
        sys.stderr.write('\nERROR: could not index BAM. Is samtools installed?\n')
