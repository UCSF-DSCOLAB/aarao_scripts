import gzip
import numpy as np
import os
import pandas as pd

# GRCh38
# wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
# mm10
# wget ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz

def get_kv_pairs(gtf_field_8):
    gf8 = gtf_field_8.rstrip(b';').split(b';')
    results = {x.split()[0].strip(): x.split()[1].strip().strip(b'"') for x in gf8 }
    return results

def ensg_to_hugo(gtf_file):
    hugo_map, pam_oguh = {}, {}
    with gzip.open(gtf_file, 'rb') as iff:
        for line in iff:
            if line.startswith(b'#'):
                continue
            line = line.strip().split(b'\t')
            if line[2] != b'gene':
                continue
            k_v = get_kv_pairs(line[8])
            hugo_map[k_v[b'gene_id'].decode("utf-8")] = k_v[b'gene_name'].decode("utf-8")
            if k_v[b'gene_name'].decode("utf-8") not in pam_oguh:
                pam_oguh[k_v[b'gene_name'].decode("utf-8")] = k_v[b'gene_id'].decode("utf-8")
    return hugo_map, pam_oguh


def csv_to_10x(csv_file, outdir, annotation, annotation_r):
    x = pd.read_csv(csv_file, index_col=0)
    x = x[sorted(x.columns)]
    with open(os.path.join(outdir, 'barcodes.tsv'), 'w') as off:
        for barcode in x.columns:
            print('{}-1'.format(barcode), file=off)
    genes = dict(zip(sorted(annotation), range(len(annotation))))
    with open(os.path.join(outdir, 'genes.tsv'), 'w') as off:
        for gene in genes:
            print('{}\t{}'.format(gene, annotation[gene]), file=off)

    with open(os.path.join(outdir, 'matrix.mtx'), 'w') as off:
        print(r'%%MatrixMarket matrix coordinate integer general',
              r'%', sep='\n', file=off)
        print(len(genes), len(x.columns), (x>0).sum().sum()+3, sep=' ', file=off)
        for tup in x.itertuples():
            try:
                gid = genes[annotation_r[tup[0]]]+1
            except KeyError:
                assert os.path.splitext(tup[0])[0] in annotation_r, tup[0]
                gid = genes[annotation_r[os.path.splitext(tup[0])[0]]]+1
            for i in range(len(tup)-1):
                if tup[i+1] == 0:
                    continue
                print(gid, i+1, tup[i+1], sep=' ', file=off)
                

def main():
    h_annot, h_annot_r = ensg_to_hugo('Homo_sapiens.GRCh38.84.gtf.gz')
    m_annot, m_annot_r = ensg_to_hugo('Mus_musculus.GRCm38.84.gtf.gz')

    for sample in os.listdir():
        if not sample.endswith('csv.gz'):
            continue
        outdir = os.path.join(os.getcwd(), os.path.splitext(os.path.splitext(sample)[0])[0])
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        csv_to_10x(sample,
                   outdir,
                   h_annot if outdir.endswith('h') else m_annot,
                   h_annot_r if outdir.endswith('h') else m_annot_r)

if __name__ == '__main__':
    main()