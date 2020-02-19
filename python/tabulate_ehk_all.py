import json
import numpy as np
import os
import pandas as pd

from collections import defaultdict

KRUMMELLAB='/Users/arjunarkalrao/krummellab'
KRUMMELLAB='/data/shared/krummellab'
ehk_threshold_folder = os.path.join(KRUMMELLAB, 'ipi/data/refs/hg38_files')
ehk_vals = {'original': json.load(open(os.path.join(ehk_threshold_folder, 'ehk_original_thresholds.json'))),
            'strict': json.load(open(os.path.join(ehk_threshold_folder, 'ehk_strict_thresholds.json')))}

ehk_genes = list(ehk_vals['original'])
results = {
    'original': defaultdict(dict),
    'strict': defaultdict(dict)
}

fix_c = {
    'liveP': 'live',
    'tcellP': 'tcell',
    'sroma': 'stroma'
}

plate_df = pd.read_csv(os.path.join(KRUMMELLAB, 'arrao/projects/new_rnaseq_pipeline/testOthers/files.tsv'), sep='\t')
for plate in plate_df.itertuples():
    df = pd.read_csv(os.path.join(plate.old_folder, 'output/plate%s/plate%s.tpm_table' % (plate.plate, plate.plate)), index_col=0, sep='\t')
    df.columns = [x[:-3].rstrip('0123456789') for x in df.columns]
    df = df.loc[ehk_genes].apply(lambda x: np.log2(x+0.1))

    for _type in results:
        _res = df.apply(lambda x: x> ehk_vals['original'][x.name], axis=1).sum().to_dict()
        for x, y in _res.items():
            if 'IPI' not in x:
                continue
            s, c = x.split('.rna.')
            c = fix_c[c] if c in fix_c else c
            if c in results[_type][s]:
                results[_type][s][c] += (',%s' % y)
            else:
                results[_type][s][c] = str(y)

pd.DataFrame.from_dict(results['original']).transpose().to_csv('all_ipi_original_ehk.tsv', sep='\t', na_rep='NA')
pd.DataFrame.from_dict(results['original']).transpose()[['live', 'myeloid', 'stroma', 'tcell', 'treg', 'tumor']].to_csv('all_ipi_original_ehk_6pops.tsv', sep='\t', na_rep='NA')
pd.DataFrame.from_dict(results['strict']).transpose().to_csv('all_ipi_strict_ehk.tsv', sep='\t', na_rep='NA')
pd.DataFrame.from_dict(results['strict']).transpose()[['live', 'myeloid', 'stroma', 'tcell', 'treg', 'tumor']].to_csv('all_ipi_strict_ehk_6pops.tsv', sep='\t', na_rep='NA')


