from __future__ import division, print_function

import argparse
import json
import logging
import numpy as np
import os
import pandas as pd
import re


from collections import defaultdict

CTRL_REGEX = re.compile(r'(?P<ctrl>(CONTROL))_'
                        r'(?P<juhr>(jurkat|UHR))[.]'
                        r'(?P<plate>(Plate))'
                        r'(?P<pltnum>([0-9]{1,2}|nova))')

#sanity
compartments = {'live',
                'tcell',
                'treg',
                'myeloid',
                'tumor',
                'stroma',
                'epcam',
                'microglia',
                'cd11bposhladrneg',
                'cd45neg',
                'cd45pos',
                'dpost',
                'tpost',
                'tpre',
                'teff',
                'thelper'
                }

IPI_REGEX = re.compile(r'^(?P<ipi>(IPI[A-Z]{2,4}[0-9]{3}))[.]'
                       r'(?P<tn>(([TNRML]|NASH)[0-9]))[.]'
                       r'rna[.]'
                       r'(?P<stain>(' + '|'.join(compartments) + '))'
                       r'(?P<puro>([Pp]*))'
                       r'(?P<suffix>([1-9]*))$')


class MyUniversalHelpFormatter(argparse.HelpFormatter):
    """
    This formatter formats both the description and argument defaults formatting
    of the argparse help string.
    """
    def _fill_text(self, text, width, indent):
        """
        This module was taken from the ArgumentDefaultsHelpFormatter class
        within argparse.  It deals with the formatting of arguments in that it
        appends the default value to teh description of each argument.
        """
        return ''.join(indent + line for line in text.splitlines(True))

    def _get_help_string(self, action):
        """
        This module was taken from the RawDescriptionHelpFormatter class
        within argparse.  It deals with the formatting of the description string
        and allows properly formatted descriptions ot be printed without line
        wrapping.
        """
        help_ = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help_ += ' (default: %(default)s)'
        return help_


def get_sample_tuple(sample):
    sample_re = IPI_REGEX.search(sample) or CTRL_REGEX.search(sample)
    assert sample_re
    if 'ctrl' in sample_re.groupdict():
      # CONTROL
      patient = '{}_{}'.format(sample_re['ctrl'], sample_re['juhr']) 
      sample = '{}_{}.{}{}'.format(sample_re['ctrl'], sample_re['juhr'], sample_re['plate'], 
                                   sample_re['pltnum'])
      indication = patient
      tissue_type = patient
      compartment = patient
    else:
      patient = sample_re['ipi']
      sample = '{}.{}'.format(sample_re['ipi'], sample_re['tn']) 
      indication = sample_re['ipi'][3:].rstrip('0123456789')
      tissue_type = sample_re['tn']
      compartment = sample_re['stain'] + sample_re['puro']
    return patient, sample, indication, tissue_type, compartment


def read_hugo_map(mapfile):
    hugo_map = pd.read_csv(mapfile, sep='\t', header=None, names=['chr', 'ensembl', 'hugo'])
    hugo_map = hugo_map.loc[~hugo_map['ensembl'].isna()]
    hugo_map.index = hugo_map['ensembl']
    hugo_map.drop('ensembl', axis=1, inplace=True)
    return hugo_map.to_dict()['hugo'], list(hugo_map[hugo_map['chr'] == 'MT'].index)


def get_fastp_stats(basedir):
    logging.info('Reading fastp metrics')
    results = pd.DataFrame(columns = [
      'raw_read_count', 'raw_base_count', 'raw_mean_length', 
      'filtered_read_count', 'filtered_base_count', 'filtered_mean_length', 
      'retention_rate', 'gc_content', 'duplication_rate', 'insert_size_peak'])
    dtypes = {'raw_read_count': int,
              'raw_base_count': int,
              # 'raw_mean_length': # Already str/object
              'filtered_read_count': int,
              'filtered_base_count': int,
              # 'filtered_mean_length': # Already str/object
              # 'retention_rate',  # Already floats
              # 'gc_content',  # Already floats
              # 'duplication_rate',  # Already floats
              'insert_size_peak': int
    }
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        with open(os.path.join(basedir, sample, 'fastp.json')) as iff:
            jsf = json.load(iff)
        results.loc[sample] = [
            jsf['summary']['before_filtering']['total_reads'],
            jsf['summary']['before_filtering']['total_bases'],
            '%s,%s' % (jsf['summary']['before_filtering']['read1_mean_length'], jsf['summary']['before_filtering']['read2_mean_length']),
            jsf['summary']['after_filtering']['total_reads'],
            jsf['summary']['after_filtering']['total_bases'],
            '%s,%s' % (jsf['summary']['after_filtering']['read1_mean_length'], jsf['summary']['after_filtering']['read2_mean_length']),
            round(100 * jsf['summary']['after_filtering']['total_bases'] / jsf['summary']['before_filtering']['total_bases'], 2),
            round(100 * jsf['summary']['after_filtering']['gc_content'], 2),
            round(100 * jsf['duplication']['rate'], 2),
            jsf['insert_size']['peak']
        ]
    for d in dtypes:
        results[d] = results[d].astype(dtypes[d])
    return results.sort_index(axis=0)


def get_rrna_reads(basedir):
    logging.info('Reading rrna counts')
    results = pd.DataFrame(columns = ['reads'])
    dtypes = {'reads': int}
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        with open(os.path.join(basedir, sample, sample + '.trimmed.rrna.sorted.flagstat')) as iff:
            results.loc[sample] = [iff.readline().split()[0]]
    for d in dtypes:
        results[d] = results[d].astype(dtypes[d])
    return results.sort_index(axis=0)


def get_star_results(basedir):
    logging.info('Reading STAR alignment metrics')
    results = pd.DataFrame(columns = ['input_reads',
                                      'uniq_map_reads', 'uniq_map_pct',
                                      'multimap_lte20_reads', 'multimap_lte20_pct',
                                      'multimapp_gt20_reads', 'multimap_gt20_pct',
                                      'unmapped_pct',
                                      'chimeric_reads','chimeric_pct'])
    dtypes = {'input_reads': int,
              'uniq_map_reads': int,
              # 'uniq_map_pct',   Already float
              'multimap_lte20_reads': int,
              # 'multimap_lte20_pct',  Already float
              'multimapp_gt20_reads': int,
              # 'multimap_gt20_pct',  Already float
              # 'unmapped_pct',  Already float
              'chimeric_reads': int,
              # 'chimeric_pct'  Already float
    }
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        out = []
        with open(os.path.join(basedir, sample,
                               sample + '.trimmed.non_rrna.star.Log.final.out')) as iff:
            lines = iff.readlines()
            # input_reads
        out.append(int(lines[5].strip().split('\t')[1]) * 2)
        # uniq_map, uniq_map_pct
        out.append(int(lines[8].strip().split('\t')[1]) * 2)
        out.append(float(lines[9].strip().split('\t')[1].rstrip('%%')))
        # multimap_lte20, multimap_lte20_pct
        out.append(int(lines[23].strip().split('\t')[1]) * 2)
        out.append(float(lines[24].strip().split('\t')[1].rstrip('%%')))
        # multimapp_gt20, multimap_gt20_pct
        out.append(int(lines[25].strip().split('\t')[1]) * 2)
        out.append(float(lines[26].strip().split('\t')[1].rstrip('%%')))
        # unmapped
        out.append(round(sum([float(l.strip().split('\t')[1][:-1]) for l in lines[28:31]]), 2))
        # chimeric, chimeric_pct
        out.append(int(lines[32].strip().split('\t')[1]) * 2)
        out.append(float(lines[33].strip().split('\t')[1].rstrip('%%')))
        results.loc[sample] = out
    for d in dtypes:
        results[d] = results[d].astype(dtypes[d])
    return results.sort_index(axis=0)


def get_duplication_percents(basedir):
    logging.info('Reading duplication metrics')
    results = pd.DataFrame(columns = ['duplication_pct'])
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        _file = os.path.join(
            basedir, sample,
            sample + '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.duplication_metrics')
        with open(_file) as iff:
            lines = iff.readlines()
        assert lines[6].startswith('LIBRARY')
        results.loc[sample] = [round(float(lines[7].strip().split()[9]) * 100, 2)]
    return results.sort_index(axis=0)


def get_rnaseq_metrics(basedir):
    logging.info('Reading rnaseq metrics')
    results = pd.DataFrame(columns = ['filter_passing_bases', 'aligned_bases',
                                      'ribo_pct', 'coding_pct', 'utr_pct', 'intron_pct', 'intergenic_pct',
                                      'median_cv_coverage', 'median_5p_bias', 'median_3p_bias'])
    dtypes = {'filter_passing_bases': int, 
              'aligned_bases': int,
              # 'ribo_pct',  Already floats
              # 'coding_pct',  Already floats
              # 'utr_pct',  Already floats
              # 'intron_pct',  Already floats
              # 'intergenic_pct', Already floats
              # 'median_cv_coverage', Already floats
              # 'median_5p_bias', Already floats
              # 'median_3p_bias' Already floats
    }
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        _file = os.path.join(
            basedir, sample,
            sample + '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated.rnaseq_metrics')
        with open(_file) as iff:
            lines = iff.readlines()
        assert lines[6].startswith('PF_BASES')
        _lines = lines[7].strip().split()
        out = []
        out.extend([int(x) for x in _lines[:2]])
        out.extend([round(float(x) * 100, 2) for x in _lines[15:20]])
        out.extend([round(float(x), 4) for x in _lines[23:26]])
        results.loc[sample] = out
    for d in dtypes:
        results[d] = results[d].astype(dtypes[d])
    return results.sort_index(axis=0)


def get_gene_expression(basedir, fpkm_genes = None):
    logging.info('Reading Gene Expression values')
    if fpkm_genes is None:
        fpkm_genes = ['ENSG00000162769',  # FLVCR1
                      'ENSG00000115486',  # GGCX
                      'ENSG00000138175',  # ARL3
                      'ENSG00000125741',  # OPA3
                      'ENSG00000155034',  # FBXL18
                      'ENSG00000196652',  # ZKSCAN5
                      'ENSG00000137413',  # TAF8
                      'ENSG00000176407',  # KCMF1
                      'ENSG00000133704']  # IPO8
    tpms = {}
    fpkms = {}
    counts = {}
    gene_exp_table = pd.DataFrame(columns=['read_counts', 'expression', 'rna_seq', 'ensembl_name'])
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        df = pd.read_csv(os.path.join(basedir, sample, sample + '.rsem.genes.results'), sep='\t',
                           index_col=0, usecols=['gene_id', 'expected_count', 'TPM', 'FPKM'])
        _df = df[['expected_count', 'TPM']].copy(deep=True)
        _df.columns = ['read_counts', 'expression']
        _df['rna_seq'] = sample
        _df['ensembl_name'] = _df.index
        gene_exp_table = gene_exp_table.append(_df, ignore_index=True, sort=False)
        tpms[sample] = df['TPM']
        counts[sample] = df['expected_count']
        fpkms[sample] = df['TPM'].loc[fpkm_genes]

    tpms = pd.concat(tpms, axis=1)    
    counts = pd.concat(counts, axis=1)    
    fpkms = pd.concat(fpkms, axis=1)    
    return (tpms.sort_index(axis=1), counts.sort_index(axis=1), fpkms.sort_index(axis=1), gene_exp_table)


def get_isoform_expression(basedir):
    logging.info('Reading Isoform Expression values')
    tpms = {}
    for sample in os.listdir(basedir):
        if not os.path.isdir(os.path.join(basedir, sample)):
            continue
        df = pd.read_csv(os.path.join(basedir, sample, sample + '.rsem.isoforms.results'), sep='\t',
                         index_col=0, usecols=['transcript_id', 'TPM'])
        tpms[sample] = df['TPM']

    tpms = pd.concat(tpms, axis=1)
    return tpms.sort_index(axis=1)


def main():
    parser = argparse.ArgumentParser(formatter_class=MyUniversalHelpFormatter)
    parser.add_argument('run_dir', help='path_to_run_dir')
    parser.add_argument('--plate', help='plate number (integer > 0 or `Nova`)', type=str, required=True)
    parser.add_argument('--transcriptome_build', help='Transcriptome build used', type=str,
                        required=False, default='Ensembl_GRCh38.85')
    parser.add_argument('--pipeline_version', help='Pipeline version', type=str,
                        required=False, default='v1.0')
    parser.add_argument('--overwrite', help='overwrite results', action='store_true')
    parser.add_argument('--logfile', type=str, default=None, help='Path to a logfile. Default is '
                        'STDOUT.')
    parser.add_argument('--log_level', type=str, choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help='The level of logging above which messages should be '
                        'printed.')
    params = parser.parse_args()

    if not (params.plate.isdigit() or params.plate == 'Nova'):
        raise RuntimeError('`--plate` can only be an integer > 0 or `Nova`')

    if params.logfile is not None:
        params.logfile = os.path.abspath(os.path.expanduser(params.logfile))
        assert os.path.exists(os.path.dirname(params.logfile)), 'Invalid Directory to write logfile'
        print('Logging to file: %s' % params.logfile)
    else:
        print('Logging to STDOUT')
    logging.basicConfig(filename=params.logfile,
                        filemode='a',
                        level=getattr(logging, params.log_level),
                        format='%(levelname)s: %(message)s')

    run_dir = os.path.abspath(params.run_dir)
    assert os.path.exists(run_dir)
    assert {'output', 'metrics', 'log'} < set(os.listdir(run_dir))
    results_dir = os.path.join(run_dir, 'results')
    if os.path.exists(results_dir):
        if not params.overwrite:
            raise RuntimeError('%s exists. Cowardly exiting without overwriting data. use '
                               '--overwrite to overwrite results.' % results_dir)
    else:
        os.mkdir(results_dir)

    fastp = get_fastp_stats(os.path.join(run_dir, 'metrics'))
    rrna = get_rrna_reads(os.path.join(run_dir, 'metrics'))
    star_results = get_star_results(os.path.join(run_dir, 'log'))
    duplication = get_duplication_percents(os.path.join(run_dir, 'metrics'))
    rnaseq_metrics = get_rnaseq_metrics(os.path.join(run_dir, 'metrics'))
    gene_tpm, gene_counts, gene_fpkm, gene_exp_table = get_gene_expression(os.path.join(run_dir,
                                                                                        'output'))
    isoform_tpm = get_isoform_expression(os.path.join(run_dir, 'output'))

    ensembl_hugo_map_file = os.path.join(os.environ['KLAB'],
                                         'ipi/data/refs/hg38_files/hugo_to_ensg_with_chr.tsv')
    hugo_map, mt_genes = read_hugo_map(ensembl_hugo_map_file)
    hugo_map = defaultdict(str, hugo_map)

    coding_counts = pd.DataFrame(
        {'chromosomal': gene_counts.loc[~gene_counts.index.isin(mt_genes)].apply(sum).astype(int), 
         'mitochondrial': gene_counts.loc[mt_genes].apply(sum).astype(int) 
         })
    coding_counts['mt_pct'] = (100 * coding_counts['mitochondrial']/(coding_counts['mitochondrial'] + coding_counts['chromosomal'])).apply(lambda x: round(x, 2))

    gene_exp_table['hugo_name'] = gene_exp_table.apply(lambda x: hugo_map[x.ensembl_name], axis=1)

    ehk_threshold_folder = os.path.join(os.environ['KLAB'], 'ipi/data/refs/hg38_files')

    ehk_files = {'EHK_legacy': os.path.join(ehk_threshold_folder, 'ehk_original_thresholds.json'),
                 'EHK_legacy_strict': os.path.join(ehk_threshold_folder, 'ehk_strict_thresholds.json'),
                 'EHK': os.path.join(ehk_threshold_folder, 'ehk_new_thresholds.json')}

    ehk_vals = {x: json.load(open(y)) for x, y in ehk_files.items()}

    ehk_all_genes = pd.read_csv(os.path.join(ehk_threshold_folder, 'EHK_all_genes.tsv'), index_col=0, header=0, sep='\t') 
    ehk_all_genes = ehk_all_genes[~ehk_all_genes['ENSG'].isna()]

    exp_qc = pd.DataFrame(index=gene_tpm.columns)

    for ev in ehk_vals:
        exp_qc[ev] = gene_tpm.apply(lambda x: sum([np.log2(x[_x]+0.1) >= ehk_vals[ev][_x]
                                                   for _x in ehk_vals[ev]]), axis=0)

    exp_qc['low_exp_flag'] = gene_fpkm.transpose().apply(
        lambda x: (sum([x[_x] == 0 for _x in ('ENSG00000162769',  # FLVCR1
                                              'ENSG00000115486',  # GGCX
                                              'ENSG00000138175',  # ARL3
                                              'ENSG00000125741',  # OPA3
                                              'ENSG00000155034')]) >= 3 or  # FBXL18
                   any(x[_x] < 1 for _x in ('ENSG00000196652',  # ZKSCAN5
                                            'ENSG00000137413',  # TAF8
                                            'ENSG00000176407',  # KCMF1
                                            'ENSG00000133704'))),  # IPO8
        axis=1)

    exp_qc['expressed_genes'] = (gene_tpm > 0.0).apply(sum, axis=0)
    exp_qc['expressed_EHK_genes'] = (gene_tpm.loc[ehk_all_genes['ENSG']] > 0.0).apply(sum, axis=0)



    info = pd.DataFrame(index=fastp.index)
    info['patient'], info['sample'], info['indication'], info['tissue_type'], info['compartment'] = \
      zip(*fastp.index.map(get_sample_tuple))
    info['sequencer'] = 'NovaSeq' if (params.plate == 'Nova' or int(params.plate) > 25) else 'HiSeq'
    info['rna_seq_plate'] = params.plate

    # Now that we have the info table, let's populate the cell counts table
    if os.path.exists(os.path.join(run_dir, 'cell_counts.tsv')):
        #cell_counts = pd.read_csv(os.path.join(run_dir, 'cell_counts.tsv'), sep='\t', header=True,
        #                          index_col=0)
        raise NotImplementedError('Can\'t process cell_counts.tsv yet')
    else:
        cell_counts = pd.DataFrame(index=info.index)
        cell_counts['well'] = 'NA'
        cell_counts['submission_name'] = 'NA'
        cell_counts['submission_date'] = 'NA'
        cell_counts['cell_count'] = 'NA'
        cell_counts['stain_panel'] = 'NA'

    assert sorted(cell_counts.index) == sorted(info.index)
    cell_counts = cell_counts.loc[info.index]

    info = pd.merge(info, cell_counts, left_index=True, right_index=True)

    info['transcriptome_build'] = params.transcriptome_build
    info['pipeline_version'] = params.pipeline_version
    

    outdict = {
        'info': info,
        'fastp': fastp,
        'ribosomal_rna': rrna,
        'star_alignment': star_results,
        'coding_counts': coding_counts,
        'MarkDuplicates': duplication,
        'CollectRnaSeqMetrics': rnaseq_metrics,
        'QC': exp_qc
    }

    logging.info('Generating RNAseq table')
    # sanity
    keys = fastp.index
    for x in outdict:
        assert all(outdict[x].index == keys)
    rnaseq_table = pd.concat(outdict, axis=1, keys=['info', 'fastp', 'ribosomal_rna',
                                                    'star_alignment', 'coding_counts',
                                                    'MarkDuplicates', 'CollectRnaSeqMetrics',
                                                    'QC'])

    # Fill in columns that require 
    rnaseq_table['QC', 'adapter_pct'] = rnaseq_table.apply(lambda x: round(100 - x['fastp', 'retention_rate'], 2), axis=1)
    rnaseq_table['QC', 'usable_map_pct'] = rnaseq_table.apply(lambda x: round(100 * x['CollectRnaSeqMetrics', 'aligned_bases']/x['fastp', 'filtered_base_count'], 2), axis=1)
    rnaseq_table['QC', 'overall_map_pct'] = rnaseq_table.apply(lambda x: round(100 * x['CollectRnaSeqMetrics', 'aligned_bases']/x['fastp', 'raw_base_count'], 2), axis=1)

    logging.info('Writing outputs')
    rnaseq_table.to_csv(os.path.join(results_dir, 'rnaseq_table.tsv'), sep='\t', header=True,
                        index=True)
    gene_tpm.to_csv(os.path.join(results_dir, 'gene_tpm_table.tsv'), sep='\t', header=True,
                    index=True)
    gene_counts.to_csv(os.path.join(results_dir, 'gene_counts_table.tsv'), sep='\t', header=True,
                       index=True)
    isoform_tpm.to_csv(os.path.join(results_dir, 'isoform_tpm_table.tsv'), sep='\t', header=True,
                       index=True)
    gene_exp_table.to_csv(os.path.join(results_dir, 'gene_exp_table.tsv'), sep='\t', header=True,
                       index=False, columns=['rna_seq', 'ensembl_name', 'hugo_name', 'read_counts',
                                             'expression'])

if __name__ == '__main__':
    main()

