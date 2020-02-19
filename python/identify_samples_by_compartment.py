import argparse
import os
import pandas as pd
import re

from collections import Counter

def main():
    parser = argparse.ArgumentParser()
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('--patients', type=str, help='A table with one patient per line. Must contain '
                     'a key `patient` in the header. All other columns will be inherited by the '
                     'output.')
    grp.add_argument('--indication', type=str, help='An entire indication to process.')
    parser.add_argument('--clinical', type=str, required=False, default=None)
    parser.add_argument('--compartments', type=str, nargs='+', required=True)
    grp2 = parser.add_mutually_exclusive_group(required=False)
    grp2.add_argument('--split_by_compartment', action='store_true', required=False)
    grp2.add_argument('--outfile', type=str, required=False, default=None)
    parser.add_argument('--ehk_cutoff', type=int, required=False, default=None)
    parser.add_argument('--rseq_table_dir', type=str,
                        default=os.path.join(os.environ['KLAB'],
                                             'arrao/data/gene_expression/rsq_tables'),
                        required=False)
    params = parser.parse_args()
    #params = parser.parse_args(['--patients', 'gyn_patients.list', '--compartments', 'live', 'myeloid', 'tcell', 'tumor', 'stroma', 'epcam'])
    #params = parser.parse_args(['--indication', 'GYN', '--compartments', 'live', 'myeloid', 'tcell', 'tumor', 'stroma', 'epcam'])
    #params = parser.parse_args(['--indication', 'GYN', '--clinical', '../histology.tsv', '--split_by_compartment', '--compartments', 'live', 'myeloid', 'tcell', 'tumor', 'stroma', 'epcam'])

    if params.patients:
        df = pd.read_csv(params.patients, sep='\t', header=0, index_col=None)
        if 'patient' not in df:
            raise RuntimeError('We need a column named patient to proceed!')

    else:
        df = pd.DataFrame({'patient': [''.join(['IPI', params.indication, '{:03d}'.format(x+1)]) 
                                       for x in range(999)]})
        params.patients = '.'

    if params.clinical:
        dfc = pd.read_csv(params.clinical, sep='\t', header=0, index_col=None)
        df = pd.merge(df, dfc, on='patient', how='left')
        df[df.isna()] = 'unknown'

    df.index = df['patient']    
    additional_groups = [g for g in df.columns if g not in ['patient']]
    if 'group' in additional_groups:
        g = [i for i, _g in enumerate(additional_groups) if _g == 'group']
        assert len(g) == 1, "Cannot have multiple columns called group"
        g = additional_groups.pop(g[0])
        additional_groups = [g] + additional_groups

    columns = ['sample', 'plate', 'compartment'] + additional_groups + ['EHK', 'low_exp_flag']
    samples = pd.DataFrame(columns=columns)

    ipi_regex = re.compile(r'^(?P<ipi>(IPI[A-Z]{2,4}[0-9]{3}))\.'
                           r'(?P<tn>(([TNRML]|NASH)[0-9]))\.'
                           r'rna\.'
                           r'(?P<compartment>([a-zA-Z]+))[2-9]*$')
    i = 0
    for rseq_table in os.listdir(params.rseq_table_dir):
        if rseq_table.startswith('.'):
            continue
        rseq_table = os.path.join(params.rseq_table_dir, rseq_table)
        if os.path.isdir(rseq_table):
            continue
        plate = re.sub('_rnaseq_table.tsv', '', os.path.basename(rseq_table))
        if plate in ('plate5', 'plateNova'):
            continue


        print('Processing %s' % rseq_table)
        rseq_table = pd.read_csv(rseq_table, sep='\t', header=[0,1], index_col=0)
        _samples = {s: ipi_regex.match(s) for s in rseq_table.index if ipi_regex.match(s)}

        for s in _samples:
            _s = _samples[s]
            if _s.group('compartment') not in params.compartments:
                continue
            if _s.group('ipi') not in df['patient']:
                continue
            out_record = [s, plate, _s.group('compartment')] + \
                         [df.loc[_s.group('ipi'), g] for g in additional_groups] + \
                         [rseq_table.loc[s, 'QC']['EHK'], rseq_table.loc[s, 'QC']['low_exp_flag']]
            samples.loc[i] = out_record
            #samples.loc[i] = [s, plate, df.loc[_s.group('ipi'), 'group'], rseq_table.loc[s, 'QC']['EHK'], rseq_table.loc[s, 'QC']['low_exp_flag']]
            i += 1

    if params.ehk_cutoff:
        samples = samples.loc[samples['EHK']>=params.ehk_cutoff]
    dupes = [x for x, y in Counter(samples['sample']).items() if y > 1]
    if dupes:
        print('WARNING: The following samples have dupes -- %s' % ', '.join(dupes))

    if params.split_by_compartment:
        for compartment in params.compartments:
            outfile = os.path.join(os.path.dirname(params.patients), '{}_valid_samples.tsv'.format(compartment))
            samples[samples['compartment']==compartment].to_csv(outfile, sep='\t', header=True, index=False)
    else:
        outfile = params.outfile or os.path.join(os.path.dirname(params.patients), 'valid_samples.tsv')
        samples.to_csv(outfile, sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()
