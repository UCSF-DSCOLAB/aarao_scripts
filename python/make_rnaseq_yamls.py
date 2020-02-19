"""
This program will parse all fastqs in a given folder and produce a yml file that is compatible with
ruby_pipe's rna_seq pipeline.

The value passed to --filename_re_order will be compiled into a regex that will be used to extract
information from the fastq file names. 
E.g. by default, --regex_oder is 
            <IPI>_<TN>_<DRNA>_<COMPARTMENT><SPACER><R12><SPACER>

and this translates into the regex

r'(?P<ipi>((IPI)?[A-Z]{2,4}[0-9]{3}))_(?P<tn>((([TN])|(NASH))[0-9]))_(?P<drna>([dr]na))_(?P<compartment>((epcam)|(live)|(myeloid)|(stroma)|(tcell)|(treg)|(tumor)|(teff)|(CD45pos))).*(?P<r12>(R[12])).*'



Usage: 
1. Naive usage
    python generate_yml.py 30

    This will create /data/shared/krummellab/ipi/input_ymls/plates/plate30.rna_paired_align.yml
    using the fastqs in /cbc2/data2/samadb/IPI/fastqs/Plate30Rnaseq and assumes that all files are
    named in the following format
           `<IPI>_<TN>_<DRNA>_<COMPARTMENT><SPACER><R12><SPACER>`

    E.g. IPIHNSC088_T1_rna_stroma_S100_L004_R2_001.fastq.gz

    The order of the regexes can be changed using --filename_re_order using the same notations as above.
    The regexes for the fields can be changed using the appropriate flags.

2. Advanced usage example
    python generate_yml.py 30 \\
                           --fastq_dir /data/fastqs/plate30 \\
                           --outfile /data/outputs/plate30 \\
                           --overwrite \\
                           --compartment_re '((tcell)|(treg)|(tumor)|(CD45_high))'

    This will attempt to create the file /data/outputs/plate30/plate30.rna_paired_align.yml if
    /data/outputs/plate30 is a directory, or will create the file /data/outputs/plate30 if it is not
    a directory. We will overwrite the file if it exists. Fastqs are pulled from
    /data/fastqs/plate30 and files are assumed to be named in the following format
           `<IPI>_<TN>_<DRNA>_<COMPARTMENT><SPACER><R12><SPACER>`

    E.g. IPIHNSC088_T1_rna_tcell_S100_L004_R2_001.fastq.gz

    `IPIHNSC088_T1_rna_stroma_S100_L004_R2_001.fastq.gz` will be flagged since there is no match for
    the compartment regex.

"""
from __future__ import print_function

import argparse
import glob
import os
import re
import yaml

def get_filename_regex(re_type, regexes, expected_tokens, template):
    """
    compile a regex for a filename given an template and a dict of regexes to replace tokens in the
    template string.

    :param str re_type: 'control' or 'fastq'.
    :param dict(str, str) regexes: A dict containing tokens expected in the remplate as keys and
           the corresponding regex as value.
    :param set expected_tokens: A set of expected tokens in the template
    :param str template: The template for the regex string.
    :return: The compiled regex
    :rtype: 

    """
    assert re_type in ['control', 'fastq']
    token_re = re.compile(r'<(?P<token>([A-Z0-9]*))>')
    tokens = {x.group('token') for x in  re.finditer(token_re, template)}
    diff = tokens - expected_tokens
    if diff:
        raise RuntimeError('Unexpected regex tokens seen in %s_re_order:  %s' %
                           (re_type, str(diff)))
    diff = expected_tokens - tokens
    if diff:
        raise RuntimeError('Expected regex tokens not seen in %s_re_order:  %s' %
                           (re_type, str(diff)))
    # Now build the regex
    out_re = template
    for token in tokens:
        out_re = re.sub('<%s>' % token, regexes[token], out_re)
    out_re = re.compile(out_re, re.IGNORECASE)
    return out_re

def process_fastqs(fastq_dir, fastq_re, control_re, plate_number):
    results = {}
    fail = False
    for fastq_file in glob.glob(os.path.join(fastq_dir,'*','*.gz')):
        _fastq_file = os.path.basename(fastq_file)
        details = re.search(fastq_re, _fastq_file) or re.search(control_re, _fastq_file)
        if not details:
            print('ERROR: %s does not match regex.' % fastq_file)
            fail = True
            continue
        ## the group either has a control, or an IPI
        if 'ipi' in details.groupdict():
            sample_name = '.'.join([details.group('ipi'),
                                    details.group('tn'),
                                    details.group('drna'),
                                    details.group('compartment')])
        else:
            sample_name = '%s.Plate%s' % (details.group('control'), plate_number)
        if sample_name not in results:
            results[sample_name] = {
                'R1': [],
                'R2': [],
            }
        results[sample_name][details.group('r12')].append(fastq_file)
    return results, fail

        

def main():
    default_fastq_dir = '/cbc2/data2/samadb/IPI/fastqs/Plate<PLATE_NUMBER>Rnaseq'
    default_outfile = '/data/shared/krummellab/ipi/input_ymls/plates/' \
                      'plate<PLATE_NUMBER>.rna_paired_align.yml'
    parser = argparse.ArgumentParser(description=str(__doc__), usage='See Below',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('plate_number', type=int, help='Plate number')
    parser.add_argument('--fastq_dir', help='Fastq Directory', default=default_fastq_dir,
                        required=False)
    parser.add_argument('--outfile', help='Output yaml', default=default_outfile,
                        required=False)
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output file?')
    parser.add_argument('--ipi_re', help='ipi regex', default=r'(IPI)?[A-Z]{2,4}[0-9]{3}',
                        required=False)
    parser.add_argument(
        '--compartment_re', help ='compartment regex.',
        default=r'(epcam)|(live)|(myeloid)|(stroma)|(tcell)|(treg)|(tumor)|(teff)',
        required=False)
    parser.add_argument('--tn_re', help ='Tumor/normal regex.', default=r'(([TN])|(NASH))[0-9]',
                        required=False)
    parser.add_argument('--control_re', help ='control regex.', default=r'CONTROL_jurkat',
                        required=False)
    parser.add_argument('--drna_re', help ='dna/rna regex.', default=r'[dr]na', required=False)
    parser.add_argument('--r12_re', help ='read 1/2 regex.', default=r'R[12]', required=False)
    parser.add_argument('--fastq_re_template', help='This is the actual order of the above ' \
                        'regexes for the fastq file name',
                        default='<IPI>_<TN>_<DRNA>_<COMPARTMENT><SPACER><R12><SPACER>',
                        required=False)
    parser.add_argument('--control_re_template', help='This is the actual order of the above ' \
                        'regexes for the control file name.',
                        default='<CONTROL><SPACER><R12><SPACER>', required=False)
    params = parser.parse_args()

    params.plate_number = str(params.plate_number)

    if params.fastq_dir == default_fastq_dir:
        params.fastq_dir = re.sub('<PLATE_NUMBER>', params.plate_number, params.fastq_dir)    
    assert os.path.exists(params.fastq_dir), 'Fastq dir %s does not exist' % params.fastq_dir


    if params.outfile == default_outfile:
        params.outfile = re.sub('<PLATE_NUMBER>', params.plate_number, params.outfile)
    else:
        if os.path.isdir(params.outfile):
            params.outfile = os.path.join(params.outfile,
                                          'plate%s.rna_paired_align.yml' % params.plate_number)

    if os.path.exists(params.outfile) and not params.overwrite:
        raise RuntimeError('Outfile exists. Please use a new name, or use --overwrite')

    regexes = {
        'IPI': "(?P<ipi>(%s))" % params.ipi_re,
        'TN': "(?P<tn>(%s))" % params.tn_re,
        'DRNA': "(?P<drna>(%s))" % params.drna_re,
        'COMPARTMENT': "(?P<compartment>(%s))" % params.compartment_re,
        'SPACER': ".*",
        'R12': "(?P<r12>(%s))" % params.r12_re,
        'CONTROL': "(?P<control>(%s))" % params.control_re
    }
    ## Process filename regex
    fastq_re = get_filename_regex(
        'fastq', regexes, {'IPI', 'DRNA', 'COMPARTMENT', 'SPACER', 'R12', 'TN'},
        params.fastq_re_template)
    contol_re = get_filename_regex(
        'control', regexes, {'CONTROL', 'SPACER', 'R12'}, params.control_re_template)


    fastq_files, fail = process_fastqs(params.fastq_dir, fastq_re, contol_re, params.plate_number)

    if fail:
        print('INFO: Please fix regexes to accommodate failed samples and try again.')
        exit(0)
    output_yaml = {
        ':cohort_name': 'plate%s' % params.plate_number,
        ':flat_reference_gtf': '/cbc2/data1/data/samadb/ipi/ref/flat_reference_GRCh38.85.gtf',
        ':frag_size': 100,
        ':genes_transcripts': './genes_transcripts.txt',
        ':genome': 'hg38',
        ':keep_temp_files': True,
        ':log_dir': './log',
        ':metrics_dir': './metrics',
        ':modules': ['rsem'],
        ':output_dir': './output',
        ':pipe': 'rna',
        ':read_size': 100,
        ':scratch_dir': './scratch',
        ':script': 'paired_align',
        ':tmp_dir': '/scratch/%s' % os.environ['USER'],
        ':verbose': True
    }
    output_samples = {':samples': []}
    for sample in fastq_files:
        assert len(fastq_files[sample]['R1']) == len(fastq_files[sample]['R2'])
        assert len(fastq_files[sample]['R1']) == 1, 'DOES NOT COMPUTE'
        _sample = {
            ':sample_name': sample,
            ':replicates': [{
                ':inputs': [{
                    ':fq1': fastq_files[sample]['R1'][0],
                    ':fq2': fastq_files[sample]['R2'][0],
                }]
            }]
        }
        output_samples[':samples'].append(_sample)
            
    with open(params.outfile, 'w') as off:
        print(yaml.dump(output_yaml, default_flow_style=False), file=off, end='')
        print(yaml.dump(output_samples, default_flow_style=False), file=off)

if __name__ == '__main__':
    main()