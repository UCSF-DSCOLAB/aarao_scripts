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


    # process cell_barcodes.json
    out_cbc_json = os.path.join('cellranger_outs', 'cell_barcodes.json')
    in_cbc_json = os.path.join(cellranger_dir, 'cell_barcodes.json')
    with open(in_cbc_json, 'r') as iff:
        _in_cbc_json = json.load(iff)
    _in_cbc_json = [c for c in _in_cbc_json if c in keep_barcodes]
    with open(out_cbc_json, 'w') as off:
        json.dump(_in_cbc_json, off, indent=4, sort_keys=True)

    shutil.copy(os.path.join(cellranger_dir, 'metrics_summary.csv'),
                os.path.join('cellranger_outs', 'metrics_summary.csv'))
    shutil.copy(os.path.join(cellranger_dir, 'web_summary.html'),
                os.path.join('cellranger_outs', 'web_summary.html'))
    if os.path.exists(os.path.join(cellranger_dir, 'vdj_reference')):
        shutil.copytree(os.path.join(cellranger_dir, 'vdj_reference'),
                        os.path.join('cellranger_outs', 'vdj_reference'))

    out_airr_tsv = os.path.join('cellranger_outs', 'airr_rearrangement.tsv')
    in_airr_tsv = os.path.join(cellranger_dir, 'airr_rearrangement.tsv')
    if os.path.exists(in_airr_tsv):
        _in_airr_tsv = pd.read_csv(in_airr_tsv, header=0, index_col=0, sep="\t")
        _in_airr_tsv = _in_airr_tsv.loc[_in_airr_tsv.index.isin(keep_barcodes)].copy()
        _in_airr_tsv.to_csv(out_clonotypes_csv, header=True, index=True, sep="\t")

    print("ACTION: Processing all_contig files ....")
    process_file_groups(cellranger_dir, keep_barcodes, prefix='all_contig', contig_str='contig', 
                        library_name=library_name, index_fasta=True, sample_num=sample_num, 
                        update_verbosity=update_verbosity)
    for f in os.listdir('all_contig_outs'):
        os.rename(os.path.join('all_contig_outs', f),
                  os.path.join('cellranger_outs', f))
    shutil.rmtree('all_contig_outs')

    print("ACTION: Processing filtered_contig files ....")
    process_file_groups(cellranger_dir, keep_barcodes, prefix='filtered_contig', 
                        contig_str='contig', library_name=library_name, index_fasta=False, 
                        sample_num=sample_num, update_verbosity=update_verbosity)
    for f in os.listdir('filtered_contig_outs'):
        os.rename(os.path.join('filtered_contig_outs', f),
                  os.path.join('cellranger_outs', f))
    shutil.rmtree('filtered_contig_outs')

    clonotypes_of_interest = pd.read_csv(os.path.join('cellranger_outs', 'all_contig_annotations.csv'), 
                                         header=0, index_col=0)
    clonotypes_of_interest = set(clonotypes_of_interest.loc[clonotypes_of_interest['raw_clonotype_id'] != "None"]["raw_clonotype_id"])

    out_clonotypes_csv = os.path.join('cellranger_outs', 'clonotypes.csv')
    in_clonotypes_csv = os.path.join(cellranger_dir, 'clonotypes.csv')
    
    _in_clonotypes_csv = pd.read_csv(in_clonotypes_csv, header=0, index_col=0)
    _in_clonotypes_csv = _in_clonotypes_csv.loc[_in_clonotypes_csv.index.isin(clonotypes_of_interest)].copy()
    _in_clonotypes_csv.to_csv(out_clonotypes_csv, header=True, index=True)
    
    print("ACTION: Processing concat_ref files ....")
    process_file_groups(cellranger_dir, keep_barcodes=keep_barcodes, prefix='concat_ref', 
                        contig_str='concat_ref', library_name=library_name, 
                        clonotypes_of_interest=clonotypes_of_interest, index_fasta=True, 
                        sample_num=sample_num, update_verbosity=update_verbosity)
    for f in os.listdir('concat_ref_outs'):
        os.rename(os.path.join('concat_ref_outs', f),
                  os.path.join('cellranger_outs', f))
    shutil.rmtree('concat_ref_outs')
    print("ACTION: Processing consensus files ....")
    process_file_groups(cellranger_dir, keep_barcodes=keep_barcodes, prefix='consensus', 
                        contig_str='consensus', library_name=library_name,
                        clonotypes_of_interest=clonotypes_of_interest, index_fasta=True, 
                        sample_num=sample_num, update_verbosity=update_verbosity)
    for f in os.listdir('consensus_outs'):
        os.rename(os.path.join('consensus_outs', f),
                  os.path.join('cellranger_outs', f))
    shutil.rmtree('consensus_outs')
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
    outdir = '{}_outs'.format(prefix)
    os.mkdir(outdir)
    # First process the fasta
    out_fasta = os.path.join(outdir, '{}.fasta'.format(prefix))
    in_fasta = os.path.join(cellranger_dir, '{}.fasta'.format(prefix))
    
    cname_re = re.compile('_{}.*'.format(contig_str))

    contigs_of_interest = clonotypes_of_interest or keep_barcodes
    contigs_to_keep = []
    with open(in_fasta, 'r') as iff:
        with open(out_fasta, 'w') as off:
            sname = None
            seq = []
            for line in iff:
                if line.startswith('>'):
                    if sname is not None:
                        CB = cname_re.sub('', sname[1:])
                        if CB in contigs_of_interest:
                            print(sname, *seq, sep='\n', file=off)
                            contigs_to_keep.append(sname[1:])
                        else:
                            # redact it in the fasta so the bam processing is easier.
                            seq = 'N' * len(''.join(seq))
                            print(sname, *seq, sep='\n', file=off)
                    sname = line.strip()
                    seq = []
                else:
                    seq.append(line.strip())
            CB = cname_re.sub('', sname[1:])
            if CB in keep_barcodes:
                print(sname, *seq, sep='\n', file=off)
                contigs_to_keep.append(sname[1:])
    
    if index_fasta:
        pysam.faidx(out_fasta)

    # Second process the annotations csv
    out_ann_csv = os.path.join(outdir, '{}_annotations.csv'.format(prefix))
    in_ann_csv = os.path.join(cellranger_dir, '{}_annotations.csv'.format(prefix))
    
    if os.path.exists(in_ann_csv):
        in_ann_df = pd.read_csv(in_ann_csv, index_col=0)
        in_ann_df = in_ann_df.loc[in_ann_df.index.isin(keep_barcodes)].copy()
        in_ann_df.to_csv(out_ann_csv, sep=',', index=True, header=True)

    # Third process the annotations bed
    out_ann_bed = os.path.join(outdir, '{}_annotations.bed'.format(prefix))
    in_ann_bed = os.path.join(cellranger_dir, '{}_annotations.bed'.format(prefix))

    if os.path.exists(in_ann_bed):
        in_ann_df = pd.read_csv(in_ann_bed, index_col=None, header=None, sep='\t')
        in_ann_df = in_ann_df.loc[in_ann_df[0].isin(contigs_to_keep)].copy()
        in_ann_df.to_csv(out_ann_bed, sep='\t', index=False, header=False)

    # Fouth process the annotations json
    out_ann_json = os.path.join(outdir, '{}_annotations.json'.format(prefix))
    in_ann_json = os.path.join(cellranger_dir, '{}_annotations.json'.format(prefix))

    if os.path.exists(in_ann_json):
        with open(in_ann_json, 'r') as iff:
            _in_ann_json = json.load(iff)
        _in_ann_json = [x for x in _in_ann_json if x['contig_name'] in contigs_to_keep]
        with open(out_ann_json, 'w') as off:
            json.dump(_in_ann_json, off, indent=4, sort_keys=True)

    # Fouth process the "annotations fastq".
    out_fastq = os.path.join(outdir, '{}.fastq'.format(prefix))
    in_fastq = os.path.join(cellranger_dir, '{}.fastq'.format(prefix))

    if os.path.exists(in_fastq):
        with open(in_fastq, 'r') as iff:
            with open(out_fastq, 'w') as off:
                # This assumes the fastq is 4 lines per record (i.e. the seq and qual strings are 
                # one-liners)
                i = 0
                contig_is_valid = False
                for line in iff:
                    if i > 0:
                        if contig_is_valid:
                            print(line, file=off, end='')
                        i -= 1
                    elif line.startswith('@'):
                        contig_is_valid = False
                        i = 3
                        if line[1:].strip() in contigs_to_keep:
                            print(line, file=off, end='')
                            contig_is_valid = True
                    else:
                        raise RuntimeError('Invalid fastq file')

    # Lastly, process the bam and input fastq. This is the file that can contain be used to generate
    # the fastqs but remember, these CBs only match teh whitelist provided.

    forward = 'ACGTN'
    reverse = 'TGCAN'
    trans = str.maketrans(forward, reverse)

    out_bam = os.path.join(outdir, '{}.bam'.format(prefix))
    in_bam = os.path.join(cellranger_dir, '{}.bam'.format(prefix))
    
    if not os.path.exists(in_bam):
        # Nothing to do here
        return None
    in_sam_handle = pysam.Samfile(in_bam, 'rb')
    out_sam_handle = pysam.Samfile(out_bam, 'wb', template=in_sam_handle)

    if prefix == 'all_contig':
        # Since we'll retain a majority of reads, it's easier to keep a list of reads we'll be 
        # dumping so we can filter the bam later
        reads_to_remove = gzip.open(os.path.join(outdir, 'reads_to_remove.list.gz'), 'wt')
    num_reads = retained_reads = 0
    try:
        for read in in_sam_handle.fetch(until_eof=True):
            num_reads += 1
            if num_reads % update_verbosity == 0:
                print('Processed {} reads'.format(num_reads))

            if prefix == 'all_contig':

                rdict = dict(read.tags)
                if 'CB' not in rdict or rdict['CB'] not in keep_barcodes:
                    if read.is_read2:
                    # protects from duplication
                        print('@{}'.format(read.qname), file=reads_to_remove)
                    continue
            else:
                if read.reference_name not in contigs_to_keep:
                    continue
            retained_reads += 1
            # Write to the output bam
            out_sam_handle.write(read)
    except:
        print('Failed on read:')
        print(read.to_dict())
        print("######")
        print(read.tags)
        raise
    finally:
        in_sam_handle.close()
        out_sam_handle.close()
        if prefix == 'all_contig':
            reads_to_remove.close()

        print('Processed:',
              'Total Reads:{}'.format(num_reads),
              'Retained reads:{}'.format(retained_reads),
              sep='\n')
    pysam.index(out_bam)
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
    assert not os.path.isdir(params.outdir)
    params.outdir = os.path.abspath(params.outdir)

    if params.tempdir == '$TMPDIR':
        params.tempdir = None
    else:
        assert os.path.exists(params.tempdir) and os.path.isdir(params.tempdir)
    tempdir = tempfile.mkdtemp(dir=params.tempdir)
    os.chdir(tempdir)

    outputs = process_10x_vdj(params.cellranger_dir,
                              params.keep_barcodes,
                              params.library_name,
                              params.sample_num,
                              params.update_verbosity)
    
    print("ACTION: Moving output to final destination....")
    shutil.move(tempdir,
                params.outdir)

if __name__ == '__main__':
    main()