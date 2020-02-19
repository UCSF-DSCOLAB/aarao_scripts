import json
import numpy as np
import os
import pandas as pd
import re

ehk_threshold_folder = os.path.join(os.environ['KLAB'], 'ipi/data/refs/hg38_files')

ehk_vals = {'original': json.load(open(os.path.join(ehk_threshold_folder, 'ehk_original_thresholds.json'))),
             'strict': json.load(open(os.path.join(ehk_threshold_folder, 'ehk_strict_thresholds.json')))}

ehk_genes = list(ehk_vals['original'].keys())

data = {}
for plate in range(26):
  plate += 1
  data[plate] = {
  'new': pd.read_csv(os.path.join(os.environ['KLAB'], 'arrao/data/gene_expression/tpms', 'plate%s_gene_tpm_table.tsv' % plate), header=0, index_col=0, sep='\t'),
  'old': pd.read_csv(os.path.join(os.environ['KLAB'], 'arrao/data/gene_expression/old_tpms', 'plate%s.tpm_table' % plate), header=0, index_col=0, sep='\t'),
  }
  for k in data[plate]:
    data[plate][k] = pd.DataFrame(data[plate][k].loc[ehk_genes]).apply(lambda x: np.log2(x+0.1))
    data[plate][k] = data[plate][k][sorted(data[plate][k])]

  data[plate]['old'].columns = [re.sub('.r0$', '', x) for x in data[plate]['old'].columns]
  
for plate in range(26):
  plate += 1
  print(plate)
  print('new_not_in_old: ', ', '.join(sorted(set(data[plate]['new'].columns) - set(data[plate]['old'].columns))))
  print('old_not_in_new: ', ', '.join(sorted(set(data[plate]['old'].columns) - set(data[plate]['new'].columns))))

