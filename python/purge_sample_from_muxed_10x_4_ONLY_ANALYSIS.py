import argparse
import glob
import gzip
import h5py
import numpy as np
import os
import pandas as pd
import pysam
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cellranger_dir', type=str, required=True, help='The path to a '
                        'cellranger-count `outs` folder.')
    parser.add_argument('--outdir', type=str, required=True, help='The path to an output directory'
                        '.')
    parser.add_argument('--keep_barcodes', type=str, required=True, help='A file containing '
                        'barcodes belonging to valid samples (whitelist) with one barcode per line '
                        'and no header.')
    params = parser.parse_args()

    print("ACTION: Validating inputs....")
    assert os.path.exists(params.cellranger_dir) and os.path.isdir(params.cellranger_dir)
    params.cellranger_dir = os.path.abspath(params.cellranger_dir)
    assert os.path.exists(params.keep_barcodes) and os.path.isfile(params.keep_barcodes)
    params.keep_barcodes = os.path.abspath(params.keep_barcodes)
    assert os.path.isdir(params.outdir)
    params.outdir = os.path.abspath(params.outdir)

    if is_gzipfile(params.keep_barcodes):
        keep_barcodes = {l.strip() for l in gzip.open(params.keep_barcodes, 'rt')}
    else:
        keep_barcodes = {l.strip() for l in open(params.keep_barcodes, 'r')}

    print("ACTION: Analysis folder....")
    os.chdir(os.path.join(params.outdir, 'cellranger_outs'))
    files_to_process = glob.glob('analysis/**/clusters.csv', recursive=True) +\
                       glob.glob('analysis/**/projection.csv', recursive=True)
    for f in files_to_process:
        _f = pd.read_csv(f, header=0, index_col=0)
        _f = _f.loc[_f.index.isin(keep_barcodes)].copy()
        _f.to_csv(f, sep=',', header=True, index=True)
    
if __name__ == '__main__':
    main()