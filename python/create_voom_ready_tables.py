import argparse
import os
import pandas as pd
import re

def ensg_to_hugo(gtf_file):
    hugo_map, pam_oguh = {}, {}
    with gzip.open(gtf_file, 'rb') as iff:
        for line in iff:
            if line.startswith(b'#'):
                continue
            line = line.strip().split(b'\t')
            if line[2] != b'gene':
                continue
            k_v = get_kv_pairs(line[8])
            hugo_map[k_v[b'gene_id'].decode("utf-8")] = k_v[b'gene_name'].decode("utf-8")
            if k_v[b'gene_name'].decode("utf-8") not in pam_oguh:
                pam_oguh[k_v[b'gene_name'].decode("utf-8")] = k_v[b'gene_id'].decode("utf-8")
    return hugo_map, pam_oguh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', type=str, required=True)
    parser.add_argument('--hugo_map', type=str, default=os.path.join(os.environ['KLAB'], 
                                                                     'ipi/data/refs/hg38_files',
                                                                     'hugo_to_ensg.tsv'))
    parser.add_argument('--counts_dir', type=str,
                        default=os.path.join(os.environ['KLAB'],
                                             'arrao/data/gene_expression/counts'),
                        required=False)
    #parser.add_argument('--blacklist', type=str, default=None, required=False)
    params = parser.parse_args()
    #params = parser.parse_args(['--samples', 'gyn_pr_myeloid.tsv'])
    samples = pd.read_csv(params.samples, sep='\t', header=0, index_col=0)
    plates = True
    if 'plate' not in samples.columns:
        samples['plate'] = ''
        plates = False

    hugo_map = pd.read_csv(params.hugo_map, sep='\t', header=0, index_col=0)
    hugo_map.drop_duplicates(keep=False, inplace=True)


    df = None
    for counts_file in os.listdir(params.counts_dir):
        plate = re.sub('_gene_counts_table.tsv', '', re.sub('[.]counts_table$', '', counts_file))
        if plate in ['plate5', 'plateNova']:
            continue
        counts_file = os.path.join(params.counts_dir, counts_file)
        if os.path.isdir(counts_file):
            continue
        print('Processing %s' % counts_file)
        if plates and plate not in samples['plate'].unique():
            continue
        _df = pd.read_csv(counts_file, sep='\t', header=0, index_col=0, nrows=0)
        _samples = [re.sub(r'\.r0', '', s) for s in _df.columns]
        if plates:
            # Assuming that if you know th eplate, you know the exact sample name
            _samples = {s: i+1 for i, s in enumerate(_samples) if s in samples.loc[samples['plate']==plate].index}
        else:
            _samples = {s: i+1 for i, s in enumerate(_samples) if s.rstrip('0123456789') in samples.index}
        if _samples:
            print('Identified %s samples to pull: %s' % (len(_samples), ', '.join(_samples)))
            for s in _samples:
                samples.loc[s, 'plate'] = plate
            _df = pd.read_csv(counts_file, sep='\t', header=0,
                              index_col=0, usecols=[0]+list(_samples.values()))
            if df is None:
                df = _df
            else:
                df = pd.merge(df, _df, left_index=True, right_index=True)
 
    if not plates:
        for sample in samples[samples.plate == ''].index:
            print('Dropping sample %s due to a lack of data' % sample)
        samples = samples[samples.plate != '']
        samples.to_csv(''.join([os.path.splitext(params.samples)[0], '_with_plate',
                                os.path.splitext(params.samples)[1]]), sep='\t', header=True,
                       index=True)

    df.columns = [re.sub('.r0$', '', s) for s in df.columns]
    df = df[samples.index]
    df.index = [(hugo_map.loc[x, 'HUGO'] if x in hugo_map.index else x) for x in df.index]
    df.to_csv(''.join([os.path.splitext(params.samples)[0], '_counts.tsv']), sep='\t', header=True,
              index=True)


if __name__ == '__main__':
    main()
