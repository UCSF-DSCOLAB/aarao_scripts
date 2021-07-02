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

def update_fastqs_and_bam(cellranger_dir, keep_barcodes, GEX_library_name, AUX_library_names, 
                          GEX_sample_num=1, AUX_sample_nums=None, update_verbosity=10000):
    """
    :param str cellranger_dir: The path to the cellranger directory with the bam and the counts 
                               matrices and h5 files
    :param str keep_barcodes: A list of barcodes meant to be reatined in the output. Assumes the 
                              barcodes are named in the same format as the bam (RE: -1 suffix)
    :param str GEX_library_name: The name of the GEX library
    :param list(str) AUX_library_names: The name(s) of any auxiliary libraries in the bam
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

    forward = 'ACGTN'
    reverse = 'TGCAN'
    trans = str.maketrans(forward, reverse)

    in_bam = os.path.join(cellranger_dir, 'possorted_genome_bam.bam')
    outputs = {'bams': ['possorted_genome_bam.bam']}

    in_sam_handle = pysam.Samfile(in_bam, 'rb')
    out_sam_handle = pysam.Samfile(outputs['bams'][0], 'wb', template=in_sam_handle)

    # in_sam_handle.header['RG'] is an ordered dict
    read_groups = {x['ID']: x['LB'] for x in in_sam_handle.header['RG']}
    outputs['fastqs'] = [GEX_library_name] + AUX_library_names
    _sample_nums = [GEX_sample_num] + _AUX_sample_nums

    libraries = {}
    sample_nums = {}
    try:
        for x in in_sam_handle.header['RG']:
            if x['LB'] in libraries:
                continue
            i = len(libraries)
            libraries[x['LB']] = outputs['fastqs'][i]
            sample_nums[x['LB']] = _sample_nums[i]
            os.mkdir(libraries[x['LB']])
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
        in_sam_handle.close()
        out_sam_handle.close()

    fastq_handles = {x: {} for x in libraries}
    
    num_reads = retained_reads = 0
    try:
        for read in in_sam_handle.fetch(until_eof=True):
            num_reads += 1
            if num_reads % update_verbosity == 0:
                print('UPDATE: Processed {} reads'.format(num_reads))

            rdict = dict(read.tags)
            if 'CB' not in rdict:
                continue
            if rdict['CB'] not in keep_barcodes:
                continue
            retained_reads += 1

            read_group = rdict['RG']
            lane = read_group[-1]
            library = read_groups[read_group]

            if lane not in fastq_handles[library]:
                fastq_handles[library][lane] = {
                    'R1': gzip.open('{}/{}_S{}_L00{}_R1_001.fastq.gz'.format(libraries[library],
                                                                             libraries[library], 
                                                                             sample_nums[library], 
                                                                             lane), 'wt'),
                    'R2': gzip.open('{}/{}_S{}_L00{}_R2_001.fastq.gz'.format(libraries[library],
                                                                             libraries[library], 
                                                                             sample_nums[library], 
                                                                             lane), 'wt'),
                    'I1': gzip.open('{}/{}_S{}_L00{}_I1_001.fastq.gz'.format(libraries[library],
                                                                             libraries[library], 
                                                                             sample_nums[library], 
                                                                             lane), 'wt')
                }

            # Write to the output bam
            out_sam_handle.write(read)

            # Write to the output fastqs
            # I1
            print('@{qname} {rnum}:N:0:{i7index}\n'
                  '{i7index}\n'
                  '+\n'
                  '{i7index_qual}'.format(qname=read.qname, 
                              rnum=1,
                              i7index=rdict['BC'],
                              i7index_qual=rdict['QT']),
                  file=fastq_handles[library][lane]['I1'])
            # R1
            print('@{qname} {rnum}:N:0:{i7index}\n'
                  '{raw_bc}{raw_umi}\n'
                  '+\n'
                  '{raw_bc_qual}{raw_umi_qual}'.format(qname=read.qname,
                                                       rnum=1,
                                                       i7index=rdict['BC'],
                                                       raw_bc=rdict['CR'],
                                                       raw_umi=rdict['UR'],
                                                       raw_bc_qual=rdict['CY'],
                                                       raw_umi_qual=rdict['UY']),
                  file=fastq_handles[library][lane]['R1'])
            # R2
            print('@{qname} {rnum}:N:0:{i7index}\n'
                  '{read_seq}\n'
                  '+\n'
                  '{read_qual}'.format(qname=read.qname, 
                              rnum=2,
                              i7index=rdict['BC'],
                              read_seq= read.seq if not read.is_reverse else read.seq.translate(trans),
                              read_qual= read.qual if not read.is_reverse else read.qual[::-1]),
                  file=fastq_handles[library][lane]['R2'])
    except:
        print('ERROR: Failed on read:')
        print(read.to_dict())
        print("######")
        print(read.tags)
        raise
    finally:
        in_sam_handle.close()
        out_sam_handle.close()
        for library in fastq_handles:
            for lane in fastq_handles[library]:
                for readname in fastq_handles[library][lane]:
                    fastq_handles[library][lane][readname].close()

        print('Update: Processed:',
              'Total Reads:{}'.format(num_reads),
              'Retained reads:{}'.format(retained_reads),
              sep='\n')
    pysam.index(outputs['bams'][0])
    outputs['bams'].append('{}.bai'.format(outputs['bams'][0]))
    return outputs

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
    #    ├─ barcode_idx                             REDACT
    #    ├─ barcode_info [HDF5 group]               
    #    │   ├─ genomes                             COPY
    #    │   └─ pass_filter                         REDACT
    #    ├─ barcodes                                COPY         
    #    ├─ count                                   REDACT
    #    ├─ feature_idx                             REDACT
    #    ├─ features [HDF5 group]                   COPY
    #    │   ├─ _all_tag_keys                     
    #    │   ├─ feature_type                      
    #    │   ├─ genome                            
    #    │   ├─ id                                
    #    │   ├─ name                              
    #    │   ├─ pattern [Feature Barcoding only]  
    #    │   ├─ read [Feature Barcoding only]     
    #    │   └─ sequence [Feature Barcoding only] 
    #    ├─ gem_group                               REDACT
    #    ├─ library_idx                             REDACT
    #    ├─ library_info                            COPY
    #    ├─ metrics [HDF5 group; see below]         COPY
    #    └─ umi                                     REDACT
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
            else if 'metrics_json' in h5_keys():
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
    print("ACTION: Processing bam....")
    outputs = update_fastqs_and_bam(params.cellranger_dir,
                                    keep_barcodes,
                                    params.GEX_library_name,
                                    params.AUX_library_names,
                                    params.GEX_sample_num,
                                    params.AUX_sample_nums,
                                    params.update_verbosity)

    os.mkdir('cellranger_outs')
    os.mkdir('fastqs')

    shutil.move(outputs['bams'][0], 'cellranger_outs')
    shutil.move(outputs['bams'][1], 'cellranger_outs')
    for library in outputs['fastqs']:
        shutil.move(library, 'fastqs/')

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

    print("ACTION: copying as-is files....")
    shutil.copytree(os.path.join(params.cellranger_dir, 'analysis'), 'cellranger_outs/analysis')
    if os.path.exists(os.path.join(params.cellranger_dir, 'feature_reference.csv')):
        shutil.copy(os.path.join(params.cellranger_dir, 'feature_reference.csv'), 'cellranger_outs/')
    shutil.copy(os.path.join(params.cellranger_dir, 'metrics_summary.csv'), 'cellranger_outs/')
    shutil.copy(os.path.join(params.cellranger_dir, 'web_summary.html'), 'cellranger_outs/')

    print("ACTION: Moving output to final destination....")
    shutil.move(os.path.join(tempdir, 'update_fastqs_and_bam'),
                params.outdir)
if __name__ == '__main__':
    main()