import argparse
import gzip
import json
import os
import pandas as pd
import pysam
import re
import tables
import tempfile
import scipy.sparse as sp_sparse
import shutil

from scipy.io import mmwrite


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

def process_10x_vdj(cellranger_dir, keep_barcodes, library_name, sample_num=1, 
                    update_verbosity=10000):
    """
    :param str cellranger_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param str keep_barcodes: A path to a file containing the barcodes meant to be reatined in the
                              output. Assumes the barcodes are named in the same format as the bam 
                              (RE: -1 suffix)
    :param str library_name: The name of the GEX library
    :param int sample_num: The SampleSheet.csv sample number for the GEX library
    :param int update_verbosity: Number of reads that must be processed for an update to be printed
    """
    os.mkdir('cellranger_outs')
    if is_gzipfile(keep_barcodes):
        keep_barcodes = {l.strip() for l in gzip.open(keep_barcodes, 'rt')}
    else:
        keep_barcodes = {l.strip() for l in open(keep_barcodes, 'r')}
    if os.path.exists(os.path.join(cellranger_dir, 'vdj_reference')):
        shutil.copytree(os.path.join(cellranger_dir, 'vdj_reference'),
                        os.path.join('cellranger_outs', 'vdj_reference'))

    out_airr_tsv = os.path.join('cellranger_outs', 'airr_rearrangement.tsv')
    in_airr_tsv = os.path.join(cellranger_dir, 'airr_rearrangement.tsv')
    if os.path.exists(in_airr_tsv):
        _in_airr_tsv = pd.read_csv(in_airr_tsv, header=0, index_col=0, sep="\t")
        _in_airr_tsv = _in_airr_tsv.loc[_in_airr_tsv.index.isin(keep_barcodes)].copy()
        _in_airr_tsv.to_csv(out_airr_tsv, header=True, index=True, sep="\t")
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cellranger_dir', type=str, required=True, help='The path to a '
                        'cellranger-count `outs` folder.')
    parser.add_argument('--outdir', type=str, required=True, help='The path to an output directory'
                        '.')
    parser.add_argument('--keep_barcodes', type=str, required=True, help='A file containing '
                        'barcodes belonging to valid samples (whitelist) with one barcode per line '
                        'and no header.')
    parser.add_argument('--library_name', type=str, required=True, help='The '
                        'name of the library (for output name) ')
    parser.add_argument('--sample_num', type=int, default=1, required=False, help='The sample '
                        'number in SampleSheet.csv to inject into the fastq name.')
    parser.add_argument('--tempdir', type=str, required=False, default="$TMPDIR", help='A temp '
                        'directory to run all operations. Must be on a device large enough to '
                        'store the bams + fastqs.')
    parser.add_argument('--update_verbosity', type=int, default=10000)
    params = parser.parse_args()

    print("ACTION: Validating inputs....")
    assert os.path.exists(params.cellranger_dir) and os.path.isdir(params.cellranger_dir)
    params.cellranger_dir = os.path.abspath(params.cellranger_dir)
    assert os.path.exists(params.keep_barcodes) and os.path.isfile(params.keep_barcodes)
    params.keep_barcodes = os.path.abspath(params.keep_barcodes)
    params.outdir = os.path.abspath(params.outdir)
    assert not os.path.exists(params.outdir)

    if params.tempdir == '$TMPDIR':
        params.tempdir = None
    else:
        assert os.path.exists(params.tempdir) and os.path.isdir(params.tempdir)
    tempdir = tempfile.mkdtemp(dir=params.tempdir)
    os.chdir(tempdir)

    process_10x_vdj(params.cellranger_dir,
                    params.keep_barcodes,
                    params.library_name,
                    params.sample_num,
                    params.update_verbosity)

    print("ACTION: Moving output to final destination....")
    shutil.copytree('cellranger_outs',
                    os.path.join(params.outdir, 'cellranger'))

    
    

if __name__ == '__main__':
    main()