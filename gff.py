from collections import defaultdict


class Feature:
    ID = 'id'
    SEP = ';'
    EQ = '='
    PARENT = 'parent'
    N = 1

    def __init__(self, line):
        self.line = line.rstrip()
        c = line.rstrip().split('\t')
        if len(c) < 9:
            raise Exception('GFF line has less than 9 columns')
        self.chro, self.origin, self.type, self.strand, attr = c[0], c[1], c[2], c[6], c[8]
        self.start, self.end, self.score, self.frame = int(c[3]), int(c[4]), c[5], c[7]
        attr = {k.lower(): v.strip('"') for k, v in
                [a.split(Feature.EQ, 1) for a in attr.split(Feature.SEP) if a]}
        self.attr = attr
        if Feature.ID in attr:
            self.name = attr[Feature.ID]
        else:
            self.name = 'noid'+Feature.N
            Feature.N += 1
        self.gene = attr['gene'] if 'gene' in attr else None
        self.parents = attr[Feature.PARENT].split(',') if Feature.PARENT in attr else []
        self.children = []

    def __eq__(self, other):
        if self is other:
            return True
        return self.chro == other.chro and self.start == other.start and self.end == other.end

    def __lt__(self, other):
        return self.start < other.start if self.is_sense() else self.end > other.end

    def __le__(self, other):
        return self.start <= other.start if self.is_sense() else self.end >= other.end

    def __gt__(self, other):
        return self.start > other.start if self.is_sense() else self.end < other.end

    def __ge__(self, other):
        return self.start >= other.start if self.is_sense() else self.end <= other.end

    def __str__(self):
        return '\t'.join([self.chro, self.origin, self.type, str(self.start), str(self.end),
                          str(self.score), self.strand, str(self.frame),
                          Feature.SEP.join([a+Feature.EQ+v for a, v in self.attr.items()])])

    def is_sense(self):
        return self.strand == '+'

    def overlaps(self, other):
        return self.chro == other.chro and (self.start <= other.start <= self.end or
                                            self.start <= other.end <= self.end or
                                            other.start <= self.start <= other.end)

    def overlaps_region(self, start, end):
        return (self.start <= start <= self.end or
                self.start <= end <= self.end or
                start <= self.start <= end)

    def distance(self, start, end):
        if self.end < start:
            return start - self.end
        if end < self.start:
            return end - self.start
        if start < self.end:
            return start - self.end
        if self.start < end:
            return self.start - end
        return start - self.end

    def contains(self, start, end):
        return self.start <= start <= end <= self.end

    def fake_parent(self):
        f = Feature(self.line)
        f.type, f.children, f.line, f.parents, f.origin = 'parent', [self.name], '', [], 'fake'
        return f


class GFF:
    filename = ""
    genes = []
    chromosomes = []
    parents = {}
    children = defaultdict(list)
    features = {}

    def __init__(self, filename, format='auto'):
        self.filename = filename
        with open(filename, 'r') as lines:
            if format == 'auto':
                if (filename.lower().endswith('.gff') or filename.lower().endswith('.gff3')):
                    format = 'gff'
                elif filename.lower().endswith('.gtf'):
                    format = 'gtf'
                else:
                    for l in lines:
                        if 'ID=' in l:
                            format = 'gff'
                            break
                        if 'gene_id ' in l:
                            format = 'gtf'
                            break
                    if format == 'auto':
                        raise Exception('Could not determine GFF/GTF format')
                    lines.seek(0)
            if format.lower() in ['gff', 'gff3']:
                Feature.SEP, Feature.EQ = ';', '='
                Feature.ID, Feature.GENE, Feature.PARENT = 'id', 'gene', 'parent'
            elif format.lower() == 'gtf':
                Feature.SEP, Feature.EQ = '; ', ' '
                Feature.ID, Feature.GENE, Feature.PARENT = 'gene_id', 'gene', 'transcript_id'
            else:
                raise Exception('Unknown Genome Annotation Format')

            last_gene = ''
            LN = 0
            for l in lines:
                LN += 1
                if l.startswith('#'):
                    continue
                f = Feature(l)
                if f.chro not in self.chromosomes:
                    self.chromosomes.append(f.chro)

                if f.type == 'gene':
                    last_gene = f.name
                    self.genes.append(f.name)
                    self.children[f.chro].append(f.name)
                else:
                    if not f.parents:
                        f.parents.append(last_gene)
                    for p in f.parents:
                        self.parents[f.name] = p
                        if p not in self.features:
                            self.features[p] = f.fake_parent()
                        self.children[p].append(f.name)
                        self.features[p].children.append(f.name)

                if f.name in self.features:
                    old = self.features[f.name]
                    f.children = old.children
                    if f.type == 'CDS':
                        f.start = min(f.start, old.start)
                        f.end = max(f.end, old.end)
                self.features[f.name] = f

    def get_genes(self, chro, strand=None, start=None, end=None):
        g = self.children[chro]
        g = [self.features[x] for x in g if strand is None or self.features[x].strand == strand]
        if start is None:
            return g
        if end is None:
            end = start
        return [x for x in g if x.start <= start <= x.end or
                x.start <= end <= x.end or
                start <= x.start <= x.end <= end]

    def get_transcripts_from_gene(self, gene):
        if isinstance(gene, str):
            gene = self.features[gene]
        return [self.features[t] for t in gene.children if self.features[t].type == 'mRNA']

    def get_exons_from_transcript(self, transcript, sort_by_direction=True):
        if isinstance(transcript, str):
            transcript = self.features[transcript]
        exons = [self.features[x] for x in transcript.children if self.features[x].type == 'exon']
        return sorted(exons, key=lambda x: x.start,
                      reverse=(not transcript.is_sense() and sort_by_direction))

    def get_unique_exons_from_gene(self, gene, sort_by_direction=True):
        if isinstance(gene, str):
            gene = self.features[gene]
        exons = {str(e.start)+' '+str(e.end): e for t in self.get_transcripts_from_gene(gene)
                 for e in self.get_exons_from_transcript(t.name)}
        return sorted(exons.values(), key=lambda x: x.start,
                      reverse=(not gene.is_sense() and sort_by_direction))

    def get_exons_from_gene(self, gene, sort_by_direction=True):
        return self.get_unique_exons_from_gene(gene, sort_by_direction)

    def get_CDS_from_transcript(self, transcript):
        return [self.features[x] for x in transcript.children if self.features[x].type == 'CDS']

    def get_annotation_from_transcript(self, transcript, pos):
        if pos < transcript.start or transcript.end < pos:
            return 'intragenic'
        if not [e for e in self.get_exons_from_transcript(transcript) if e.start <= pos <= e.end]:
            return 'intronic'
        cds = self.get_CDS_from_transcript(transcript)
        if not cds:
            return 'exonic'
        start = min([c.start for c in cds])
        end = max([c.end for c in cds])
        if pos < start and transcript.is_sense() or not transcript.is_sense() and end < pos:
            return '5UTR'
        if pos < start and not transcript.is_sense() or transcript.is_sense() and end < pos:
            return '3UTR'
        return 'CDS'

    def get_genes_by_feature(self, feature):
        return self.get_genes(feature.chro, feature.start, feature.end, feature.strand)

    def overlap(self, start, end, genes):
        return [self.features[g] for g in genes if start <= self.features[g].start < end or
                start < self.features[g].end <= end or
                self.features[g].start <= start <= end <= self.features[g].end]

    def is_exonic(self, chro, start, end=None, strand=None):
        if end is None:
            end = start
        for g in self.overlap(start, end, self.get_genes(chro, strand)):
            for t in g.children:
                for e in t.children:
                    if e.type == 'exon' and e.start <= start <= e.end:
                        start_found = True
                    if e.type == 'exon' and e.start <= end <= e.end:
                        end_found = True
        return start_found & end_found

    def get_transcripts(self, chro, start, end=None, strand=None):
        if end is None:
            end = start
        # return [f for f in self.features.values() if f.type == 'mRNA' and
        # f.start <= start <= f.end or f.start <= end <= f.end]
        return [self.features[t] for g in self.get_genes(chro, strand, start, end)
                for t in self.features[g.name].children
                if self.features[t].type in ['mRNA', 'transript']]

# return only exons of one transcript if start and end are part of the transcripts exons
    def get_exons(self, chro, start, end, strand=None):
        best = []
        best_t = None
        for g in self.get_genes(chro, strand, start, end):
            for t in self.features[g.name].children:
                ret, start_exon, end_exon = [], False, False
                for en in self.features[t].children:
                    e = self.features[en]
                    if e.type == 'exon' and (start <= e.start <= end or
                                             start <= e.end <= end):
                        ret.append(e)
                        if e.start <= start <= e.end:
                            start_exon = True
                        if e.start <= end <= e.end:
                            end_exon = True

                if start_exon & end_exon:
                    best = sorted(ret, key=lambda x: x.start)
                    best_t = self.features[t]
                    if start == best[0].start and end == best[-1].end:
                        return (best, self.features[t])
        return (best, best_t)

    def exon_start(self, chro, pos):
        return [e for e in self.get_exons(chro, pos, pos)
                if (e.sense and e.start == pos) or (not e.sense and e.end == pos)]

    def exon_end(self, chro, pos):
        return [e for e in self.get_exons(chro, pos, pos)
                if (not e.sense and e.start == pos) or (e.sense and e.end == pos)]

    def get_left_intron(self, feature):
        last = 0
        exons = [self.features[e] for p in feature.parents for e in self.features[p].children]
        exons = [e for e in exons if e.type == 'exon']
        for e in sorted(exons, key=lambda x: x.start):
            if e.name == feature.name:
                if last == 0:
                    return None
                return feature._replace(type='intron', start=last, end=e.start-1)
            last = e.end + 1

    def get_right_intron(self, feature):
        last = 0
        exons = [self.features[e] for p in feature.parents for e in self.features[p].children]
        exons = [e for e in exons if e.type == 'exon']
        for e in sorted(exons, key=lambda x: x.start, reverse=True):
            if e.name == feature.name:
                if last == 0:
                    return None
                return feature._replace(type='intron', start=e.end+1, end=last)
            last = e.start - 1

    def get_closest_feature(self, chro, start, end=None, strand=None):
        if chro not in self.chromosomes:
            return None
        if end is None:
            end = start
        s, e, offset = start, end, 0
        while(True):
            s = s - offset if offset < s else 1
            e += offset
            genes = self.get_genes(chro, strand, s, e)
            if genes:
                if len(genes) == 1:  # only feature
                    return genes[0]
                genes = sorted(genes, key=lambda x: x.distance(start, end))
                dist = genes[0].distance(start, end)
                if dist >= 0:  # closest feature (if multiple at same distance, pick any)
                    return genes[0]
                genes = [g for g in genes if g.distance(start, end) == dist]
                if len(genes) == 1:  # feature with most overlap
                    return genes[0]
                return sorted(genes, key=lambda x: x.end-x.start)[0]  # give the shortest feature
            offset += 1000  # next round with wider window
