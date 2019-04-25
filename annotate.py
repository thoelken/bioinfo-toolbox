#!/usr/bin/env python3
from gff import GFF
import sys

sys.stderr.write('parsing GFF...')
genome = GFF(sys.argv[1])
sys.stderr.write('DONE!\n')
sys.stderr.flush()
for line in sys.stdin:
    line = line.rstrip()
    c = line.split('\t')
    chro, start, end = c[0], int(c[1]), int(c[2])
    strand = None if len(c) < 6 else c[5]
    gene = genome.get_closest_feature(chro, start, end, strand)
    if gene is None:
        print(line+'\t.\t.\t0\t.\t.')
        continue
    annot_start, annot_end = 'genic', 'genic'
    for t in genome.get_transcripts_from_gene(gene):
        annot_start = genome.get_annotation_from_transcript(t, start)
        annot_end = genome.get_annotation_from_transcript(t, end)
    dist, sym = gene.distance(start, end), gene.attr['gene_symbol']
    print('%s\t%s\t%s\t%d\t%s\t%s' % (line, gene.name, sym, dist, annot_start, annot_end))
