from __future__ import print_function
import os
import re
from collections import defaultdict

def tokenize(x):
    x = re.sub('"', '', x).split('; ')
    return defaultdict(str, [_x.strip(' ;').split() for _x in x])


hg38_folder = os.path.join(os.environ['KLAB'], 'ipi/data/refs/hg38_files/')

with open(os.path.join(hg38_folder, 'Homo_sapiens.GRCh38.85.gtf')) as iff:
    with open(os.path.join(hg38_folder, 'hugo_to_ensg_with_chr.tsv'), 'w') as off1:
        with open(os.path.join(hg38_folder, 'hugo_to_ensg.tsv'), 'w') as off2:
            print('CHROM', 'ENSEMBL', 'HUGO', sep='\t', file=off1)
            print('ENSEMBL', 'HUGO', sep='\t', file=off2)
            for line in iff:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    tokens = tokenize(line[8])
                    print(line[0], tokens['gene_id'], tokens['gene_name'], sep='\t', file=off1)
                    print(tokens['gene_id'], tokens['gene_name'], sep='\t', file=off2)
                    