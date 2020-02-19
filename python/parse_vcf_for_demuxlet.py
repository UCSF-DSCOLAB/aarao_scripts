import argparse
import gzip
import os
import re

from collections import Counter
from file_methods import is_gzipfile
from gtf_methods import read_annotation_from_gtf_by_chrom
from vcf_methods import VCFRecord

plink_chrom_mapper = {
    '0': None,
    '1': '1',
    '2': '2',
    '3': '3',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '7',
    '8': '8',
    '9': '9',
    '10': '10',
    '11': '11',
    '12': '12',
    '13': '13',
    '14': '14',
    '15': '15',
    '16': '16',
    '17': '17',
    '18': '18',
    '19': '19',
    '20': '20',
    '21': '21',
    '22': '22',
    '23': 'X',
    '24': 'Y',
    '25': None,
    '26': 'MT',
}


def validate_input_genomes(gtf, vcf, genome):
    gzip_gtf = is_gzipfile(gtf) 
    open_func = gzip.open if gzip_gtf else open
    
    with open_func(gtf) as iff:
        for line in iff:
            if line.startswith('#'):
                if line.startswith('#!genome-version'):
                    if genome != line.strip().split()[1]:
                        raise RuntimeError('Annotations file genome {} does not match provided '
                            'genome {}'.format(line.strip().split()[1], genome))
                    else:
                        print('INFO: Genome in gtf file matches reported genome.')
                    break

    gzip_vcf = is_gzipfile(vcf) 
    print('INFO: Plink generated vcfs don\'t have genome in their header.')
    return gzip_gtf, gzip_vcf


def get_gene(vcf_record, annotations):
    h = len(annotations)
    l = 0
    m = h // 2

    while True:
        if annotations[m].start <= vcf_record.POS:
            if vcf_record.POS <= annotations[m].end:
                # Converged on solution
                return annotations[m].gene_name
            else:
                l = m
                m = (l + h) // 2
        else:
            h = m
            m = (l + h) // 2
        if l == m  or l == h:
            print('WARNING: Not found: {}.'.format(v))
            break
    return None


def parse_vcf(vcf_file, out_file, gzip_vcf, annotations):
    open_func = gzip.open if gzip_vcf else open
    first = True
    with open_func(vcf_file) as iff:
        with open_func(out_file, 'w' + ('b' * gzip_vcf)) as off:
            for line in iff:
                line = line.strip()
                if line.startswith('#'):
                    if line.startswith('##FORMAT') and first:
                        first = False
                        print('##FORMAT=<ID=GN,Number=1,Type=String,Description="Gene affected">', 
                              file=off)
                    print(line, file=off)
                else:
                    v = VCFRecord(line)
                    v.CHROM = plink_chrom_mapper[v.CHROM]
                    if v.CHROM is None:
                        continue
                    if len(set(v.GENOTYPES)) == 1:
                        # All the same genotype so skip
                        continue
                    else:
                        c = Counter(v.GENOTYPES)
                        if c['1/1'] + c['0/1'] > 1:
                            # More than one sample have the same genotype
                            continue
                        else:
                            g = get_gene(v, annotations[v.CHROM])
                            if g is None:
                                continue
                            else:
                                v.INFO += ';GN={}'.format(g)
                                print(v, file=off)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', help='The vcf file produced by Plink')
    parser.add_argument('gtf', help='The gtf file used to create 10X indexes')
    parser.add_argument('--genome', choices=['GRCh37', 'GRCh38', 'GRCm38'], help='Genome being '
                        'processed', default='GRCh38')
    parser.add_argument('--outfile', help='File to write to', default='<VCF>_parsed.vcf')
    parser.add_argument('--overwrite', help='Overwrite existing vcf?', action='store_true')
    #params = parser.parse_args()
    params = parser.parse_args(['/Users/arjunarkalrao/projects/XGI1/first_pool/genotyping/LH_10X/20200114_OEE64to65_KeepKattahLab_converted.vcf', '/Users/arjunarkalrao/important_data/10x/gtfs/hg38/genes.gtf'])

    assert os.path.exists(params.vcf)
    assert os.path.exists(params.gtf)

    gzip_gtf, gzip_vcf = validate_input_genomes(params.gtf, params.vcf, params.genome)

    if params.outfile == '<VCF>_parsed.vcf':
        extension_re = r'(?P<extension>(([.]vcf){0,1}([.]gz){0,1})$)'
        extension = re.search(extension_re, params.vcf)
        if extension is not None:
            params.outfile = re.sub(extension['extension']+'$', 
                                    '_parsed'+extension['extension'], 
                                    params.vcf)
        else:
            params.outfile = params.vcf + '_parsed.vcf' + ('.gz' * gzip_vcf)

    if os.path.exists(params.outfile):
        if params.overwrite:
            print('INFO: Overwriting existing output file: {}'.format(params.outfile))
        else:
            raise RuntimeError('Output file {} exists. Use --overwrite to  overwrite it else use '
                               '--outfile to specify a different name.'.format(params.outfile))
    
    open_func = gzip.open if gzip_gtf else open
    with open_func(params.gtf) as iff:
        annotations = read_annotation_from_gtf_by_chrom(iff, 'gene')

    annotations = {x: sorted(y.values(), key=lambda z: [z.start, z.end])
                   for x, y in annotations.items()
                   if x in plink_chrom_mapper.values()}
    
    parse_vcf(params.vcf, params.outfile, gzip_gtf, annotations)
