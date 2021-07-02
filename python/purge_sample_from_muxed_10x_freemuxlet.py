import argparse
import gzip
import os
import pandas as pd
import re
import tempfile
import shutil

from itertools import chain


def is_gzipfile(filename):
    """
    Attempt to ascertain the gzip status of a file based on the "magic signatures" of the file.
    This was taken from the stack overflow post
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress
    :param str filename: A path to a file
    :return: True if the file appears to be gzipped else false
    :rtype: bool
    """
    assert os.path.exists(filename), 'Input {} does not point to a file.'.format(filename)
    with open(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == b'\x1f\x8b\x08':
            return True
        else:
            return False


def process_clust1_samples_gz(freemuxlet_dir, keep_clusters, library_name):
    """
    :param str freemuxlet_dir: The path to the freemuxlet directory
    :param set(str)|list(str) keep_clusters: A set or list of clusters meant to be retained by this 
                                             script.
    :param str library_name: The name of the GEX library
    """
    os.mkdir('clust1_samples_gz_outs')
    
    in_files = os.listdir(freemuxlet_dir)

    print("ACTION: Processing .clust1.samples.gz ....")
    out_clust1_samples_gz = os.path.join('clust1_samples_gz_outs', 
                                         '{}.clust1.samples.gz'.format(library_name))
    
    in_clust1_samples_gz = [x for x in in_files if x.endswith('.clust1.samples.gz')]
    assert len(in_clust1_samples_gz) == 1
    in_clust1_samples_gz = os.path.join(freemuxlet_dir, in_clust1_samples_gz[0])


    _in_clust1_samples_gz = pd.read_csv(in_clust1_samples_gz, sep="\t", header=0, index_col=None)
    _in_clust1_samples_gz = \
        _in_clust1_samples_gz.loc[_in_clust1_samples_gz['BEST.GUESS'].apply(
                lambda x: set(x.split(',')).issubset(keep_clusters))].copy()
    
    _in_clust1_samples_gz['INT_ID'] = range(_in_clust1_samples_gz.shape[0])
    _in_clust1_samples_gz.to_csv(out_clust1_samples_gz, sep="\t", header=True, index=False)
    return set(_in_clust1_samples_gz['BARCODE'])

def process_freemuxlet(freemuxlet_dir, keep_barcodes, library_name, skip_clust1_samples_gz=True,
                       keep_clusters=None):
    """
    :param str freemuxlet_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param list(str)|set(str) keep_barcodes: A list of barcodes meant to be retained by this script.
    :param str library_name: The name of the GEX library
    :param bool skip_clust1_samples_gz: Should we skip processing .clust1.samples.gz
    :param list(str)|set(str) keep_clusters: A list of clusters meant to be retained by this script.
    """
    os.mkdir('freemuxlet_outs')

    in_files = os.listdir(freemuxlet_dir)

    # First process the cel.gz file so we can get droplet IDS to keep
    print("ACTION: Processing .cel.gz ....")
    out_cel_gz = os.path.join('freemuxlet_outs', '{}.cel.gz'.format(library_name))
    in_cel_gz = [x for x in in_files if x.endswith('.cel.gz')]
    assert len(in_cel_gz) == 1
    in_cel_gz = os.path.join(freemuxlet_dir, in_cel_gz[0])
    _in_cel_gz = pd.read_csv(in_cel_gz, sep="\t", header=0, index_col=None)
    _in_cel_gz = _in_cel_gz.loc[_in_cel_gz['BARCODE'].isin(keep_barcodes)].copy()
    keep_droplet_ids = {x: i for i, x in enumerate(_in_cel_gz['#DROPLET_ID'])}
    _in_cel_gz['#DROPLET_ID'] = range(_in_cel_gz.shape[0])
    _in_cel_gz.to_csv(out_cel_gz, sep="\t", header=True, index=False)

    # clust1.samples.gz if required (and if present)
    out_clust1_samples_gz = os.path.join('freemuxlet_outs',
                                         '{}.clust1.samples.gz'.format(library_name))
    in_clust1_samples_gz = [x for x in in_files if x.endswith('.clust1.samples.gz')]
    assert len(in_clust1_samples_gz) <= 1
    if len(in_clust1_samples_gz) == 1:
        in_clust1_samples_gz = os.path.join(freemuxlet_dir, in_clust1_samples_gz[0])
        if not skip_clust1_samples_gz:
            print("ACTION: Processing .clust1.samples.gz ....")
            _in_clust1_samples_gz = pd.read_csv(in_clust1_samples_gz, sep="\t", header=0, 
                                                index_col=None)
            _in_clust1_samples_gz = \
                _in_clust1_samples_gz.loc[_in_clust1_samples_gz['BARCODE'].isin(keep_barcodes)].copy()
            _in_clust1_samples_gz['INT_ID'] = range(_in_clust1_samples_gz.shape[0])
            _in_clust1_samples_gz.to_csv(out_clust1_samples_gz, sep="\t", header=True, index=False)
        if keep_clusters is None:
            _in_clust1_samples_gz = pd.read_csv(in_clust1_samples_gz, sep="\t", header=0, 
                                                index_col=None)
            keep_clusters = [x.split(',') for x in _in_clust1_samples_gz['BEST.GUESS'].unique()]
            keep_clusters = set(chain(*keep_clusters))

    # .lmix if present
    out_lmix = os.path.join('freemuxlet_outs', '{}.lmix'.format(library_name))
    in_lmix = [x for x in in_files if x.endswith('.lmix')]
    assert len(in_lmix) <= 1
    if len(in_lmix) == 1:
        in_lmix = os.path.join(freemuxlet_dir, in_lmix[0])
        print("ACTION: Processing .lmix ....")
        _in_lmix = pd.read_csv(in_lmix, sep="\t", header=0, index_col=None)
        _in_lmix = _in_lmix.loc[_in_lmix['BARCODE'].isin(keep_barcodes)].copy()
        _in_lmix['INT_ID'] = range(_in_lmix.shape[0])
        _in_lmix.to_csv(out_lmix, sep="\t", header=True, index=False)

    # .umi.gz
    out_umi_gz = os.path.join('freemuxlet_outs', '{}.umi.gz'.format(library_name))
    in_umi_gz = [x for x in in_files if x.endswith('.umi.gz')]
    assert len(in_umi_gz) <= 1
    if len(in_umi_gz) == 1:
        in_umi_gz = os.path.join(freemuxlet_dir, in_umi_gz[0])
        print("ACTION: Processing .umi.gz ....")
        with gzip.open(in_umi_gz, 'rb') as iff:
            with gzip.open(out_umi_gz, 'wb') as off:
                for line in iff:
                    line = line.split(b'\t')
                    if int(line[0]) in keep_droplet_ids:
                        line[0] = str.encode(str(keep_droplet_ids[int(line[0])]))
                        off.write(b'\t'.join(line))

    # .plp.gz
    print("ACTION: Processing .plp.gz ....")
    out_plp_gz = os.path.join('freemuxlet_outs', '{}.plp.gz'.format(library_name))
    in_plp_gz = [x for x in in_files if x.endswith('.plp.gz')]
    assert len(in_plp_gz) == 1
    in_plp_gz = os.path.join(freemuxlet_dir, in_plp_gz[0])
    _in_plp_gz = pd.read_csv(in_plp_gz, sep="\t", header=0, index_col=None, 
                             dtype={'#DROPLET_ID': int,
                                    'SNP_ID': int,
                                    'ALLELES': str,
                                    'BASEQS': str})
    _in_plp_gz = _in_plp_gz.loc[_in_plp_gz['#DROPLET_ID'].isin(keep_droplet_ids)].copy()
    _in_plp_gz['#DROPLET_ID'] = _in_plp_gz['#DROPLET_ID'].apply(lambda x: keep_droplet_ids[x])
    keep_snp_ids = set(_in_plp_gz['SNP_ID'])
    _in_plp_gz.to_csv(out_plp_gz, sep="\t", header=True, index=False)

    # .var.gz is untouched as it only contains info on the input regions. However, we need to parse
    # it to identify snps in vcf.gz to retain
    print("ACTION: Processing .var.gz ....")
    in_var_gz = [x for x in in_files if x.endswith('.var.gz')]
    assert len(in_var_gz) == 1

    shutil.copy(os.path.join(freemuxlet_dir, in_var_gz[0]),
                os.path.join('freemuxlet_outs', '{}.var.gz'.format(library_name)))
    var_gz = pd.read_csv(os.path.join('freemuxlet_outs', '{}.var.gz'.format(library_name)),
                         header=0, index_col=None, sep="\t")
    var_gz = var_gz.loc[var_gz['#SNP_ID'].isin(keep_snp_ids)].copy()
    var_gz = var_gz[['CHROM', 'POS', 'REF', 'ALT']].copy()
    var_gz.columns = ['vgzCHROM', 'vgzPOS', 'vgzREF', 'vgzALT']

    # clust1.vcf.gz
    out_clust1_vcf_gz = os.path.join('freemuxlet_outs', '{}.clust1.vcf.gz'.format(library_name))
    in_clust1_vcf_gz = [x for x in in_files if x.endswith('.clust1.vcf.gz')]
    assert len(in_clust1_vcf_gz) <= 1
    if len(in_clust1_vcf_gz) == 1:
        in_clust1_vcf_gz = os.path.join(freemuxlet_dir, in_clust1_vcf_gz[0])
        _in_clust1_vcf_gz = pd.read_csv(in_clust1_vcf_gz, header=None, index_col=None, sep="\t",
                                        comment="#")
        with gzip.open(in_clust1_vcf_gz, 'rt') as iff:
            with gzip.open(out_clust1_vcf_gz, 'wt') as off:
                last_comment_line = None
                for line in iff:
                    if line.startswith("#"):
                        if last_comment_line is not None:
                            # Always print the previous line read. This way the last comment line 
                            # (i.e.) the header for teh vcf won't get written here.
                            print(last_comment_line, file=off, end='')
                        last_comment_line = line
                    else:
                        break
                _in_clust1_vcf_gz.columns = last_comment_line.strip().split('\t')
                keep_columns = [x for x in _in_clust1_vcf_gz.columns 
                                if not x.startswith('CLUST')]
                keep_columns.extend(['CLUST{}'.format(x) for x in sorted(keep_clusters, 
                                                                         key=lambda x: int(x))])
                final_clust1_vcf_gz = pd.merge(_in_clust1_vcf_gz, var_gz,
                                               left_on=['#CHROM', 
                                                        'POS', 
                                                        'REF', 
                                                        'ALT'],
                                               right_on=['vgzCHROM', 
                                                         'vgzPOS', 
                                                         'vgzREF', 
                                                         'vgzALT'],
                                                         how='inner')
                assert len(final_clust1_vcf_gz.index) == len(var_gz.index)
                final_clust1_vcf_gz = final_clust1_vcf_gz[keep_columns].copy()
                final_clust1_vcf_gz.to_csv(off, sep="\t", header=True, index=False)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freemuxlet_dir', type=str, required=True, help='The path to a folder '
                        'containing freemulxetresults on a single sample.')
    parser.add_argument('--outdir', type=str, required=True, help='The path to an output directory'
                        '.')
    bc_cl = parser.add_mutually_exclusive_group(required=False)
    bc_cl.add_argument('--keep_clusters', type=str, default=None, help='A file containing '
                       'clusters meant to be retained by this script (whitelist) with '
                       'one cluster per line and no header.')
    bc_cl.add_argument('--keep_barcodes', type=str, default=None, help='A file containing '
                       'barcodes meant to be retained by this script (whitelist) with '
                       'one barcode per line and no header.')
    parser.add_argument('--library_name', type=str, required=True, help='The '
                        'name of the library (also the prefix for all files) ')
    parser.add_argument('--tempdir', type=str, required=False, default="$TMPDIR", help='A temp '
                        'directory to run all operations. Must be on a device large enough to '
                        'store the bams + fastqs.')
    params = parser.parse_args()

    print("ACTION: Validating inputs....")
    assert os.path.exists(params.freemuxlet_dir) and os.path.isdir(params.freemuxlet_dir)
    params.freemuxlet_dir = os.path.abspath(params.freemuxlet_dir)
    if params.keep_barcodes is not None:
        assert os.path.exists(params.keep_barcodes) and os.path.isfile(params.keep_barcodes)
        params.keep_barcodes = os.path.abspath(params.keep_barcodes)
    else:
        assert os.path.exists(params.keep_clusters) and os.path.isfile(params.keep_clusters)
        params.keep_clusters = os.path.abspath(params.keep_clusters)
    assert not os.path.isdir(params.outdir)
    params.outdir = os.path.abspath(params.outdir)

    if params.tempdir == '$TMPDIR':
        params.tempdir = None
    else:
        assert os.path.exists(params.tempdir) and os.path.isdir(params.tempdir)
    tempdir = tempfile.mkdtemp(dir=params.tempdir)
    os.chdir(tempdir)


    if params.keep_barcodes is not None:
        if is_gzipfile(params.keep_barcodes):
            keep_barcodes = {l.strip() for l in gzip.open(params.keep_barcodes, 'rt')}
        else:
            keep_barcodes = {l.strip() for l in open(params.keep_barcodes, 'r')}
        keep_clusters = None
    else: 
        keep_clusters = {str(x.strip()) for x in open(params.keep_clusters)}
        keep_barcodes = process_clust1_samples_gz(params.freemuxlet_dir,
                                                  keep_clusters,
                                                  params.library_name)


    process_freemuxlet(params.freemuxlet_dir,
                       keep_barcodes,
                       params.library_name,
                       skip_clust1_samples_gz=params.keep_clusters is None,
                       keep_clusters=keep_clusters)
    
    if params.keep_barcodes is None:
        os.rename(os.path.join('clust1_samples_gz_outs', 
                               '{}.clust1.samples.gz'.format(params.library_name)),
                  os.path.join('freemuxlet_outs', 
                               '{}.clust1.samples.gz'.format(params.library_name)))
        shutil.rmtree('clust1_samples_gz_outs')

    
    print("ACTION: Moving output to final destination....")
    shutil.move(os.path.join(tempdir, 'freemuxlet_outs'),
                params.outdir)

if __name__ == '__main__':
    main()