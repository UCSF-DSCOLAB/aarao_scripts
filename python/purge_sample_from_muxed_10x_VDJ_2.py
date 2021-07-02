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
    if is_gzipfile(keep_barcodes):
        keep_barcodes = {l.strip() for l in gzip.open(keep_barcodes, 'rt')}
    else:
        keep_barcodes = {l.strip() for l in open(keep_barcodes, 'r')}

    print("ACTION: Processing all_contig files ....")
    process_file_groups(cellranger_dir, keep_barcodes, prefix='all_contig', contig_str='contig', 
                        library_name=library_name, index_fasta=True, sample_num=sample_num, 
                        update_verbosity=update_verbosity)
    return None

def process_file_groups(cellranger_dir, keep_barcodes, prefix, contig_str, library_name,
                        clonotypes_of_interest=None, index_fasta=True, sample_num=1, 
                        update_verbosity=10000):
    """
    :param str cellranger_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param list(str) keep_barcodes: A list of barcodes meant to be reatined in the output. Assumes 
                                      the barcodes are named in the same format as the bam 
                                      (RE: -1 suffix)
    :param str prefix: The prefix for the files (all_contig, consensus, etc)
    :param str contig_str: The string used to refer to each contig in the files
    :param str library_name: The name of the GEX library
    :param list(str)|None clonotypes_of_interest: A list of clonotypes meant to be reatined in the 
                                                  output. Assumes the barcodes are named in the same 
                                                  format as the bam (RE: -1 suffix)
    :param bool index_fasta: Should the fasta be indexed?
    :param int sample_num: The SampleSheet.csv sample number for the GEX library
    :param int update_verbosity: Number of reads that must be processed for an update to be printed
    """
    # Lastly, process the bam and input fastq. This is the file that can contain be used to generate
    # the fastqs but remember, these CBs only match teh whitelist provided.
    in_bam = os.path.join(cellranger_dir, '{}.bam'.format(prefix))
    
    assert os.path.exists(in_bam)

    in_sam_handle = pysam.Samfile(in_bam, 'rb')

    # Since we'll retain a majority of reads, it's easier to keep a list of reads we'll be 
    # dumping so we can filter the bam later
    reads_to_remove = gzip.open('reads_to_remove.list.gz', 'wt')
    num_reads = retained_reads = 0
    try:
        for read in in_sam_handle.fetch(until_eof=True):
            num_reads += 1
            if num_reads % update_verbosity == 0:
                print('Processed {} reads'.format(num_reads))

            rdict = dict(read.tags)
            if 'CB' not in rdict or rdict['CB'] not in keep_barcodes:
                if read.is_read2:
                # protects from duplication
                    print('@{}'.format(read.qname), file=reads_to_remove)
                continue

            retained_reads += 1
    except:
        print('Failed on read:')
        print(read.to_dict())
        print("######")
        print(read.tags)
        raise
    finally:
        in_sam_handle.close()
        reads_to_remove.close()

        print('Processed:',
              'Total Reads:{}'.format(num_reads),
              'Retained reads:{}'.format(retained_reads),
              sep='\n')
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
    shutil.copy('reads_to_remove.list.gz',
                os.path.join(params.outdir, 'cellranger_outs', 'reads_to_remove.list.gz'))

    
    

if __name__ == '__main__':
    main()