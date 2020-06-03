import itertools 

class GFFRecord(object):
    # Attributes of any GFF record
    # https://en.wikipedia.org/wiki/General_feature_format#GFF_general_structure
    GFFFATTRIBUTES = [('sequence', str),
                      ('source', str),
                      ('feature', str),
                      ('start', int),
                      ('end', int),
                      ('score', str),
                      ('strand', str),
                      ('phase', str),
                      ('attributes', str)]

    def __init__(self, line, prefix_chr=False):
        """
        Converts a GFF record into a python object.

        :param str|list line: One line from a GFF file (string, or string split by tabs)
        :param bool prefix_chr: Prefix 'chr' to the sequence name?
        """
        if isinstance(line, str):
            line = line.strip().split('\t')
        if len(line) != 9:
            raise ValueError('Malformed GFF line. Must have the following '
                             'attributes: {}'.format('\n'.join(self.GFFFATTRIBUTES)))
        for (attr, func), value in zip(self.GFFFATTRIBUTES, line):
            setattr(self, attr, func(value))
        for feature in self.attributes.split(';'):
            try:
                attr, value = feature.strip().split()
                if attr in self.__dict__:
                    if isinstance(getattr(self, attr), list):
                        setattr(self, attr, getattr(self, attr) + [value.replace('"', '')])
                    else:
                        setattr(self, attr, [getattr(self, attr), value.replace('"', '')])
                else:
                    setattr(self, attr, value.replace('"', ''))
            except ValueError:
                pass
        if prefix_chr:
            self.sequence = 'chr' + self.sequence
        self.record_length = self.end - self.start + 1


    def __repr__(self):
        msg = ('GFF(sequence:{sequence}, '
               'gene:{gene_name}, '
               'feature:{feature}, ')

        if getattr(self, 'transcript_id', None) is not None:
            msg += 'transcript_id:{transcript_id}, '

        if getattr(self, 'exon_number', None) is not None:
            msg += 'exon_number:{exon_number}, '

        msg += 'start:{start}, length:{record_length})'
        return msg.format(**self.__dict__)

    def __hash__(self):
        if self.feature == 'gene':
            return hash(self.gene_id)
        elif self.feature == 'transcript':
            return hash(self.transcript_id)
        elif self.feature == 'CDS':
            return hash('CDS' + self.exon_id)
        elif self.feature == 'exon':
            return hash(self.exon_id)

    def __eq__(self, other):
        try:
            if self.feature == 'gene':
                return self.gene_id == other.gene_id
            elif self.feature == 'transcript':
                return self.transcript_id == other.transcript_id
            elif self.feature == 'CDS':
                    return other.feature == 'CDS' and (self.exon_id == other.exon_id)
            elif self.feature == 'exon':
                return self.exon_id == other.exon_id
            else:
                return False
        except AttributeError:
            return False

    def get(self, value, default=None):
        if value in self.__dict__:
            return self.__dict__[value]
        else:
            return default


def read_annotation_from_gtf(gtf_handle, annotation='gene', key='gene_name', prefix_chr=False):
    """
    Read the gene annotations into a dict
    :param file gtf_handle: A file handle to the annotation file.
    :param str annotation: An annotation type to pull from the gtf
    :param str key: The GTF feature to use as the key for the returned dictionary
    :param bool prefix_chr: Prefix 'chr' to the sequence name?
    :returns:  A dict of a gtf record for each record of type `annotation`
    :rtype: dict(string, GFFRecord)
    """
    annotations = {}
    for line in gtf_handle:
        if line.startswith('#'):
            continue
        else:
            gtf = GFFRecord(line, prefix_chr=prefix_chr)
            if gtf.feature == annotation:
                annotations[gtf.get(key)] = gtf
    return annotations

def read_annotation_from_gtf_by_chrom(gtf_handle, annotation='gene', prefix_chr=False):
    """
    Read the gene annotations into a dict of dicts
    :param file gtf_handle: A file handle to the annotation file.
    :param str annotation: An annotation type to pull from the gtf
    :param bool prefix_chr: Prefix 'chr' to the sequence name?
    :returns:  A dict of a gtf record for each record of type `annotation`
    :rtype: dict(string, GFFRecord)
    """
    annotations = {}
    for line in gtf_handle:
        if line.startswith('#'):
            continue
        else:
            gtf = GFFRecord(line, prefix_chr=prefix_chr)
            if gtf.feature == annotation:
                if gtf.sequence in annotations:
                    annotations[gtf.sequence][gtf.gene_name] = gtf
                else:
                    annotations[gtf.sequence] = {gtf.gene_name: gtf}
    return annotations


def interval_overlap(interval1, interval2):
    """
    Do 2 intervals overlap?

    :param list(chr|int, int, int) interval1: A list of (chr, start, stop) for interval 1
    :param list(chr|int, int, int) interval2: A list of (chr, start, stop) for interval 2

    :returns: True or False for the question "Does interval 2 overlap with interval1"
    :rtype: bool
    """
    if interval1[0] != interval2[0]:
        # If the interval is in a different chrom there can't be an overlap
        overlap = False
    elif interval2[1] > (interval1[2] + 1):
        # Intervals are on the same chrom. If there is a gap between the intervals, there is no
        # overlap. GTFs are 1-based and inclusive so two intervals one bp apart are considered
        # an "overlap" since they can be chained.
        overlap = False
    else:
        # Intervals are on the same chrom and they overlap.
        overlap = True
    return overlap

def collapse_intervals(intervals):
    """
    Collapse intervals into overlapping ones

    :param list(chr|int, int, int) intervals: A list of (chr, start, stop) per interval

    :returns:  A list of (chr, start, stop) per collapsed interval
    :rtype: list(tuple(str|int, int, int))
    """
    results = []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1], x[2]))
    curr_record = None
    for interval in intervals:
        if curr_record is None:
            # Base case
            curr_record = interval
        elif interval_overlap(curr_record, interval):
            # Update the "end" coordinate and go to the next interval.
            # If the new interval is a subset of the existing one, we should handle that.
            curr_record[2] = max(interval[2], curr_record[2])
        else:
            results.append(curr_record)
            curr_record = interval
    # Write the last record
    results.append(curr_record)
    return results


def get_introns(gtf_handle, prefix_chr=False, collapse_across_genes=True, only_protein_coding=True):
    """
    Get a list of intronic coordinates in the input gtf

    The way to do this is to first get overlapping exons per gene and the taking the diff. Collapsing
    across genes will correct for genes on the antisense strand that are overlapping with genes on
    the sense one.

    :param file gtf_handle: A file handle to the annotation file.
    :param bool prefix_chr: Prefix 'chr' to the sequence name?
    :param bool collapse_across_genes: Collapse introns to only keep the smallest overlaps across
                                       genes?
    :returns:  A list of intronic locations in the
    :rtype: list(tuple(str|int, int, int))
    """
    gene_records = []
    for gene, records in get_whole_gene_records(gtf_handle, annotations=['gene', 'exon']):
        if only_protein_coding:
            _gene = [r for r in records if r.feature == 'gene']
            assert(len(_gene)) == 1
            if _gene[0].gene_biotype != 'protein_coding':
                continue
        records = [[r.sequence, r.start,  r.end] for r in records if r.feature == 'exon']
        gene_record = {'exons': collapse_intervals(records)}
        # Remember, exons is ordered by the start position
        gene_record['span'] = [gene_record['exons'][0][0],  # Chrom can be from either
                               gene_record['exons'][0][1],  # Start of the first exon
                               gene_record['exons'][-1][2]]  # End of the last exon
        gene_records.append(gene_record)

    if collapse_across_genes:
        gene_records2 = sorted(gene_records, key=lambda x: x['span'])
        gene_records = []
        curr_record = None
        for gr in gene_records2:
            if curr_record is None:
                curr_record = [gr]
            elif interval_overlap(curr_record[-1]['span'], gr['span']):
                curr_record.append(gr)
            else:
                intervals = list(itertools.chain(*[cr['exons'] for cr in curr_record]))
                intervals = collapse_intervals(intervals)
                gene_records.append({'exons': intervals})  # The gene name is useless now
                curr_record = [gr]
        # Write the last record
        intervals = list(itertools.chain(*[cr['exons'] for cr in curr_record]))
        intervals = collapse_intervals(intervals)
        gene_records.append({'exons': intervals})

    introns = []
    for gr in gene_records:
        prev_exon = None
        for exon in gr['exons']:
            if prev_exon is None:
                pass
            else:
                introns.append([exon[0],  # chrom can be from either
                                prev_exon[2]+1,  # 1 base after the previous exon end
                                exon[1]-1  # one base less than the current exon start
                                ])
            prev_exon = exon

    return introns


def get_whole_gene_records(gtf_handle, annotations=['gene', 'exon'], key='gene_id',
                           prefix_chr=False):
    """
    Return one record at a time

    :param file gtf_handle: A file handle to the annotations file.
    :param str|list annotations: Annotation types to pull from the gtf
    :param str key: The GTF feature to key records by. Must be gene_name or gene_id
    :param bool prefix_chr: Prefix 'chr' to the sequence name?
    :yields: A tuple of a list of gtf record for each record of type in `annotations`
    :ytype: tuple(string, list(GFFRecord))
    :returns:  A tuple of list of a gtf record for each record of type in `annotations`
    :rtype: tuple(string, list(GFFRecord))
    """
    if isinstance(annotations, str):
        annotations = [annotations]
    elif isinstance(annotations, list):
        pass
    else:
        raise RuntimeError('annotations must be list or string')

    prev_key = None
    results = None
    for line in gtf_handle:
        if line.startswith('#'):
            continue
        else:
            gtf = GFFRecord(line, prefix_chr=prefix_chr)
            if gtf.feature in annotations:
                if gtf.get(key) == prev_key:
                    results.append(gtf)
                else:
                    if results is not None:
                        yield prev_key, results
                    prev_key = gtf.get(key)
                    results = [gtf]
    yield prev_key, results
