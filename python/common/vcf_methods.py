
complement = str.maketrans('ACTGN', 'TGACN')

class VCFRecord(object):
    # Attributes of any VCF4.0 record
    # https://www.internationalgenome.org/wiki/Analysis/vcf4.0/
    VCFATTRIBUTES = [('CHROM', str),
                     ('POS', int),
                     ('ID', str),
                     ('REF', str),
                     ('ALT', str),
                     ('QUAL', str),
                     ('FILTER', str),
                     ('INFO', str)]

    def __init__(self, line, prefix_chr=False):
        """
        Converts a VCF record into a python object
        :param str|list line: One line from a VCF file (string, or string split 
                              by tabs)
        :param bool prefix_chr: Prefix 'chr' to the sequence name?
        """
        if isinstance(line, str):
            line = line.strip().split('\t')
        if len(line) < 8:
            raise ValueError('Malformed VCF line. Must have at least the '
                             'following attributes: {}'.format(
                                    '\n'.join(self.VCFATTRIBUTES)))
        for (attr, func), value in zip(self.VCFATTRIBUTES, line[:8]):
            setattr(self, attr, func(value))
        if len(line) > 8:
            if len(line) < 9:
                raise ValueError('If a VCF record is greater than 8 fields, '
                                 'field 9 is FORMAT and there must be at least '
                                 '1 genotype field (total=10). Received '
                                 '{} fields : {}'.format(len(line), 
                                                         '|'.join(line)))
            self.FORMAT = line[8]
            self.GENOTYPES = line[9:]
    
    def __repr__(self):
        msg = ('{' + 
               '}\t{'.join([v for v,_ in self.VCFATTRIBUTES]) + 
               '}').format(**self.__dict__)
        if getattr(self, 'FORMAT', None) is not None:
            msg += '\t{}\t{}'.format(self.FORMAT, '\t'.join(self.GENOTYPES))
        return msg

    def __eq__(self, other):
        try:
            if (self.CHROM, self.POS) != (other.CHROM, other.POS):
                return False
            if self.REF == other.REF:
                return self.ALT == other.ALT
            elif self.REF == other.REF.translate(complement):
                return self.ALT == other.ALT.translate(complement)
            else:
                return False
        except AttributeError:
            return False
