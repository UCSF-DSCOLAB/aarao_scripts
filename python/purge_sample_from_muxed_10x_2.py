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
              "      --nthreads ${PBS_NUM_PPN} \\",
              "      possorted_genome_bam.bam \\",
              "      ./fastqs",
              "",
              "# Move outputs to final location",
              "mv /scratch/${{USER}}/{OF} {outdir}/{POF}".format(
                    OF=outfile[:-3], outdir=outdir, POF=outfile[8:-3]),
              sep="\n",
              file=off
              )
    return outfile

def update_h5_matrix(cellranger_dir, matrix_type, keep_barcodes):
    """
    Creates a file in the working directory with the same name as `cellranger_h5` but with the 
    values for all barcodes not in `keep_barcodes` set to 0.

    :param str cellranger_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param str matrix_type: raw or filtered?
    :param str keep_barcodes: A path to a file containing the barcodes meant to be retained in the
                              output. Assumes the barcodes are named in the same format as the bam 
                              (RE: -1 suffix)
    """
    assert matrix_type in ('raw', 'filtered')

    h5_filename = '{}_feature_bc_matrix.h5'.format(matrix_type)
    cellranger_h5 = os.path.join(cellranger_dir, h5_filename)
    shutil.copy(cellranger_h5, h5_filename)

    mtx_filename = os.path.join('{}_feature_bc_matrix'.format(matrix_type), 'matrix.mtx.gz')
    cellranger_mtx = os.path.join(cellranger_dir, mtx_filename)
    
    os.mkdir(os.path.dirname(mtx_filename))
    shutil.copy(os.path.join(os.path.dirname(cellranger_mtx), 'features.tsv.gz'),
                os.path.join(os.path.dirname(mtx_filename), 'features.tsv.gz'))
    shutil.copy(os.path.join(os.path.dirname(cellranger_mtx), 'barcodes.tsv.gz'),
                os.path.join(os.path.dirname(mtx_filename), 'barcodes.tsv.gz'))

    with gzip.open(cellranger_mtx, 'rt') as iff:
        iff.readline()
        file_comment = iff.readline().strip()[1:]  # removes the leading %

    with tables.open_file(h5_filename, 'r+') as f:
        try:
            group = f.get_node(f.root, 'matrix')
        except tables.NoSuchNodeError:
            print("Matrix group does not exist in this file.")
            return None
        barcodes = getattr(group, 'barcodes').read()
        
        for i, barcode in enumerate(barcodes):
            if barcode.decode('utf-8') not in keep_barcodes:
                group.data[group.indptr[i]:group.indptr[i+1]] = 0
        group.data.flush()

        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        with gzip.open(mtx_filename, 'wb') as off:
            mmwrite(target=off,
                    a=matrix,
                    comment=file_comment,
                    field='integer')

    return h5_filename, os.path.dirname(mtx_filename)

def update_molecule_info_h5(cellranger_dir, keep_barcodes):
    """
    Creates a molecule_info.h5 in the working directory without any barcodes not in `keep_barcodes`.

    :param str cellranger_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param str keep_barcodes: A path to a file containing the barcodes meant to be retained in the
                              output. Assumes the barcodes are named in the same format as the h5 
                              (RE: -1 suffix)
    """
    #    (root)
    #    |
    #    ├─ barcode_idx                             (266122701,) DONE
    #    ├─ barcode_info [HDF5 group]               
    #    │   ├─ genomes                             UNTOUCHED
    #    │   └─ pass_filter                         NEEDS_REDACT
    #    ├─ barcodes                                UNTOUCHED         
    #    ├─ count                                   (266122701,) DONE
    #    ├─ feature_idx                             (266122701,) DONE
    #    ├─ features [HDF5 group]                   UNTOUCHED
    #    │   ├─ _all_tag_keys                     
    #    │   ├─ feature_type                      
    #    │   ├─ genome                            
    #    │   ├─ id                                
    #    │   ├─ name                              
    #    │   ├─ pattern [Feature Barcoding only]  
    #    │   ├─ read [Feature Barcoding only]     
    #    │   └─ sequence [Feature Barcoding only] 
    #    ├─ gem_group                               (266122701,) DONE
    #    ├─ library_idx                             (266122701,) DONE
    #    ├─ library_info                            UNTOUCHED
    #    ├─ metrics [HDF5 group; see below]         UNTOUCHED
    #    └─ umi                                     (266122701,) DONE
    h5_filename = 'molecule_info.h5'
    cellranger_h5 = os.path.join(cellranger_dir, h5_filename)
    shutil.copy(cellranger_h5, '{}.bk'.format(h5_filename))

    # USEFUL: https://stackoverflow.com/a/23388605/3542991
    with h5py.File('{}.bk'.format(h5_filename), 'r') as h5_in:
        with h5py.File(h5_filename, 'w') as h5_out:
            h5_in.copy('barcodes', h5_out['/'])
            h5_in.copy('features', h5_out['/'])
            h5_in.copy('library_info', h5_out['/'])
            if 'metrics' in h5_in.keys():
                h5_in.copy('metrics', h5_out['/'])
            elif 'metrics_json' in h5_in.keys():
                h5_in.copy('metrics_json', h5_out['/'])
            else:
                assert False

            keep_bc_idxs = [i for i, b in enumerate(h5_in['barcodes']) if b.decode('utf-8') in keep_barcodes]
            mask = np.isin(h5_in['barcode_idx'], keep_bc_idxs)

            temp = h5_in['barcode_idx'][:]
            temp = temp[mask]
            h5_out.create_dataset("barcode_idx", data=temp)
            
            temp = h5_in['count'][:]
            temp = temp[mask]
            h5_out.create_dataset("count", data=temp)
            
            temp = h5_in['feature_idx'][:]
            temp = temp[mask]
            h5_out.create_dataset("feature_idx", data=temp)
            
            temp = h5_in['gem_group'][:]
            temp = temp[mask]
            h5_out.create_dataset("gem_group", data=temp)
            
            temp = h5_in['library_idx'][:]
            temp = temp[mask]
            h5_out.create_dataset("library_idx", data=temp)
            
            temp = h5_in['umi'][:]
            temp = temp[mask]
            h5_out.create_dataset("umi", data=temp)

            h5_out.create_group('barcode_info')
            h5_in.copy('barcode_info/genomes', h5_out['/barcode_info'])
            
            temp = h5_in['barcode_info']['pass_filter'][:,:]
            temp = temp[np.isin(temp[:,0], keep_bc_idxs),:]
            h5_out.create_dataset("barcode_info/pass_filter", data=temp)
    os.remove('{}.bk'.format(h5_filename))
    return h5_filename


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

    os.mkdir('cellranger_outs')
    shutil.move(out_script, 'cellranger_outs')


    print("ACTION: Processing raw h5....")
    h5_outputs = update_h5_matrix(params.cellranger_dir,
                                  'raw',
                                  keep_barcodes)
    shutil.move(h5_outputs[0], 'cellranger_outs/')
    shutil.move(h5_outputs[1], 'cellranger_outs/')

    print("ACTION: Processing filtered h5....")
    h5_outputs = update_h5_matrix(params.cellranger_dir,
                                  'filtered',
                                  keep_barcodes)
    shutil.move(h5_outputs[0], 'cellranger_outs')
    shutil.move(h5_outputs[1], 'cellranger_outs')

    print("ACTION: Processing molecule_info h5....")
    h5_outputs = update_molecule_info_h5(params.cellranger_dir, 
                                         # molecule_info h5 don't use the gem group
                                         keep_barcodes_no_gem_group)  
    shutil.move(h5_outputs, 'cellranger_outs')

    print("ACTION: Analysis folder....")
    shutil.copytree(os.path.join(params.cellranger_dir, 'analysis'), 'cellranger_outs/analysis')
    files_to_process = glob.glob('analysis/**/clusters.csv', recursive=True) +\
                       glob.glob('analysis/**/projection.csv', recursive=True)
    for f in files_to_process:
        _f = pd.read_csv(f, header=0, index_col=0)
        _f = _f.loc[_f.index.isin(keep_barcodes)].copy()
        _f.to_csv(f, sep=',', header=True, index=True)
    
    print("ACTION: copying as-is files....")
    if os.path.exists(os.path.join(params.cellranger_dir, 'feature_reference.csv')):
        shutil.copy(os.path.join(params.cellranger_dir, 'feature_reference.csv'), 'cellranger_outs/')
    shutil.copy(os.path.join(params.cellranger_dir, 'metrics_summary.csv'), 'cellranger_outs/')
    shutil.copy(os.path.join(params.cellranger_dir, 'web_summary.html'), 'cellranger_outs/')

    print("ACTION: Moving output to final destination....")
    shutil.move(os.path.join(tempdir, 'update_fastqs_and_bam'),
                params.outdir)
if __name__ == '__main__':
    main()