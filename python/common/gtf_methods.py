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
        Converts a GFF record into a python object
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
        msg = 'GFF(sequence:{sequence}, gene:{gene_name}, feature:{feature}, '

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

def read_annotation_from_gtf(gtf_handle, annotation='gene'):
    """
    Read the gene annotations into a dict
    :param file gtf_handle: A file handle ot the annotation file.
    :param str annotation: An annotation type to pull from the gtf
    :returns:  A dict of a gtf record for each record of type `annotation`
    :rtype: dict(string, GFFRecord)
    """
    annotations = {}
    for line in gtf_handle:
        if line.startswith('#'):
            continue
        else:
            gtf = GFFRecord(line)
            if gtf.feature == annotation:
                annotations[gtf.gene_name] = gtf
    return annotations

def read_annotation_from_gtf_by_chrom(gtf_handle, annotation='gene'):
    """
    Read the gene annotations into a dict of dicts
    :param file gtf_handle: A file handle ot the annotation file.
    :param str annotation: An annotation type to pull from the gtf
    :returns:  A dict of a gtf record for each record of type `annotation`
    :rtype: dict(string, GFFRecord)
    """
    annotations = {}
    for line in gtf_handle:
        if line.startswith('#'):
            continue
        else:
            gtf = GFFRecord(line)
            if gtf.feature == annotation:
                if gtf.sequence in annotations:
                    annotations[gtf.sequence][gtf.gene_name] = gtf
                else:
                    annotations[gtf.sequence] = {gtf.gene_name: gtf}
    return annotations
