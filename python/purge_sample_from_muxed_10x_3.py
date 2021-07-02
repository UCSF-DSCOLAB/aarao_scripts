import argparse
import gzip
import h5py
import numpy as np
import os
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

def update_fastqs_and_bam(cellranger_dir, keep_barcodes_file, GEX_library_name, AUX_library_names,
                          outdir, GEX_sample_num=1, AUX_sample_nums=None):
    """
    Create a bash script in the same directory to process the bams and fastqs with external tools

    :param str cellranger_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param str keep_barcodes_file: A file containing a list of barcodes meant to be reatined in the 
                                   output. Assumes the barcodes are named in the same format as the 
                                   bam (RE: -1 suffix)
    :param str GEX_library_name: The name of the GEX library
    :param list(str) AUX_library_names: The name(s) of any auxiliary libraries in the bam
    :param str outdir: A path to the output directory
    :param int GEX_sample_num: The SampleSheet.csv sample number for the GEX library
    :param list(int) AUX_sample_nums: The SampleSheet.csv sample number for all auxiliary libraries
    :param int update_verbosity: Number of reads that must be processed for an update to be printed
    """
    assert isinstance(AUX_library_names, list), "AUX_library_names must be a list (can be empty)"
    if AUX_sample_nums is None:
        _AUX_sample_nums = [GEX_sample_num + i + 1 for i in range(len(AUX_library_names))]
    else:
        assert isinstance(AUX_sample_nums, list), "AUX_sample_nums must be a list (can be empty)"
        _AUX_sample_nums = AUX_sample_nums

    in_bam = os.path.join(cellranger_dir, 'possorted_genome_bam.bam')

    in_sam_handle = pysam.Samfile(in_bam, 'rb')

    # in_sam_handle.header['RG'] is an ordered dict
    read_groups = {x['ID']: x['LB'] for x in in_sam_handle.header['RG']}
    outputs = [GEX_library_name] + AUX_library_names
    _sample_nums = [GEX_sample_num] + _AUX_sample_nums

    libraries = {}
    sample_nums = {}
    try:
        for x in in_sam_handle.header['RG']:
            if x['LB'] in libraries:
                continue
            i = len(libraries)
            libraries[x['LB']] = outputs[i]
            sample_nums[x['LB']] = _sample_nums[i]
    except IndexError:
        raise RuntimeError('Issue mapping library names and sample numbers to ({}) detected '
                           'libraries in the bam. Received GEX_library_name: ({}), '
                           'AUX_library_names: ({}), GEX_sample_num: ({}), and AUX_sample_nums: '
                           '({})'.format(len(set(read_groups.values())), 
                                         GEX_library_name,
                                         ','.join(AUX_library_names),
                                         GEX_sample_num,
                                         ','.join([str(x) for x in AUX_sample_nums]) 
                                            if AUX_sample_nums is not None else ''))
    finally:
        in_sam_handle.close()
    outfile = 'process_{}_bam_fastqs.sh'.format(GEX_library_name)
    with open(outfile, 'w') as off:
        print("#!/usr/bin/env bash",
              "#PBS -V",
              "#PBS -N {}".format(outfile[:-3]),
              "#PBS -l nodes=1:ppn=12",
              "#PBS -l vmem=150G",
              "",
              "source /krummellab/data1/ipi/software/samtools/samtools-1.10-usr/SOURCE_THIS",
              "set -e",
              "set -o nounset",
              "",
              "# Go to a temp working directory",
              "mkdir /scratch/${{USER}}/{OF}".format(OF=outfile[:-3]),
              "cd /scratch/${{USER}}/{OF}".format(OF=outfile[:-3]),
              "",
              "# Write the expected library info to a file for future reference",
              "echo \"LB : sample_name : sample_num\" > LB_info.txt",
              sep="\n",
              file=off)
        for lib in libraries:
            print("echo \"{} : {} : {}\" >> LB_info.txt".format(lib, libraries[lib], sample_nums[lib]),
                  sep="\n",
                  file=off
                  )
        print("",
              "# Extract barcodes of interest with samtools",
              "samtools view -@ ${PBS_NUM_PPN} \\",
              "              -D CB:<(zcat {}) \\".format(keep_barcodes_file),
              "              -o possorted_genome_bam.bam \\",
              "              {}".format(in_bam),
              "",
              "# Index the bamfile",
              "samtools index -@ ${PBS_NUM_PPN} possorted_genome_bam.bam",
              "",
              "# Extract fastqs from the bam with the cellranger rust tool",
              "/krummellab/data1/ipi/software/cellranger_bamtofastq/v1.3.2/bamtofastq_linux \\",
              "--nthreads ${PBS_NUM_PPN} \\",
              "possorted_genome_bam.bam \\",
              "./fastqs",
              "",
              "# Move outputs to final location",
              "mv /scratch/${{USER}}/{OF} {outdir}/{POF}".format(
                    OF=outfile[:-3], outdir=outdir, POF=outfile[8:-3]),
              sep="\n",
              file=off
              )
    return outfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cellranger_dir', type=str, required=True, help='The path to a '
                        'cellranger-count `outs` folder.')
    parser.add_argument('--outdir', type=str, required=True, help='The path to an output directory'
                        '.')
    parser.add_argument('--keep_barcodes', type=str, required=True, help='A file containing '
                        'barcodes belonging to valid samples (whitelist) with one barcode per line '
                        'and no header.')
    parser.add_argument('--GEX_library_name', type=str, required=True, help='The '
                        'name of the gene expression library (for output name) ')
    parser.add_argument('--AUX_library_names', nargs='+', type=str, default=None, required=False, 
                        help='The names of all auxiliary assays in order of appearance in the read '
                        'groups of the bam (for output folder name). Use the following command if '
                        'in doubt: `samtools view -H possorted_genome_bam.bam | tail -5`')
    parser.add_argument('--GEX_sample_num', type=int, default=1, required=False, help='The sample '
                        'number in SampleSheet.csv to inject into the GEX fastq name.')
    parser.add_argument('--AUX_sample_nums', nargs='+', type=int, default=['GEX_sample_num+n'], 
                        required=False, help='The sample numbers in SampleSheet.csv to inject into '
                        'all auxiliary fastq names.')
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
    assert not os.path.isdir(params.outdir)
    params.outdir = os.path.abspath(params.outdir)

    if params.tempdir == '$TMPDIR':
        params.tempdir = None
    else:
        assert os.path.exists(params.tempdir) and os.path.isdir(params.tempdir)
    tempdir = tempfile.mkdtemp(dir=params.tempdir)
    os.mkdir(os.path.join(tempdir, 'update_fastqs_and_bam'))
    os.chdir(os.path.join(tempdir, 'update_fastqs_and_bam'))

    if params.AUX_sample_nums == ['GEX_sample_num+n']:
        params.AUX_sample_nums = None  # update_fastqs_and_bam will handle this
    if params.AUX_library_names is None:
        params.AUX_library_names = []

    if is_gzipfile(params.keep_barcodes):
        keep_barcodes = {l.strip() for l in gzip.open(params.keep_barcodes, 'rt')}
    else:
        keep_barcodes = {l.strip() for l in open(params.keep_barcodes, 'r')}

    keep_barcodes_no_gem_group = {x.split('-')[0] for x in keep_barcodes}
    print("ACTION: Generating bash script to process bam....")
    out_script = update_fastqs_and_bam(params.cellranger_dir,
                                       keep_barcodes_file=os.path.abspath(params.keep_barcodes),
                                       GEX_library_name=params.GEX_library_name,
                                       AUX_library_names=params.AUX_library_names,
                                       outdir=params.outdir,
                                       GEX_sample_num=params.GEX_sample_num,
                                       AUX_sample_nums=params.AUX_sample_nums)

    shutil.move(out_script, os.path.dirname(params.outdir))
if __name__ == '__main__':
    main()