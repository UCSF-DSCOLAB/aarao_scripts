import glob
import os
import pandas as pd
import re

from collections import defaultdict

def rename_script(plate, old_name, new_name):
    return (
       'python /data/shared/krummellab/arrao/scripts/python/rename_rnaseq_sample.py '
       '--plate_dir /data/shared/krummellab/ipi/sequencing/processed/RSQ/{plate}_rnaseq_new/ '
       '--plate {plate} '
       '--old_sample_name {old_name} '
       '--new_sample_name {new_name}'.format(plate=plate, old_name=old_name, new_name=new_name))


CTRL_REGEX = re.compile(r'CONTROL_(jurkat|UHR).Plate([0-9]{1,2}|Nova)')
BAD_CTRL_REGEX = re.compile(r'(?P<ctrl>(control))_'
                            r'(?P<juhr>(jurkat|uhr))[.]'
                            r'(?P<plate>(plate))'
                            r'(?P<pltnum>([0-9]{1,2}|nova))', re.IGNORECASE)
good_compartments = {'live',
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

bad_compartments = {'total': 'live',
                    'CD45pos': 'cd45pos',
                    'CD45neg': 'cd45neg',
                    'cd45': 'cd45pos',
                    'Tpre': 'tpre',
                    'Tpost': 'tpost',
                    'TPo': 'tpost',
                    'Dpre': 'dpre',
                    'Dpost': 'dpost',
                    'sroma': 'stroma',
                    'livea': 'live',
                    'liveb': 'live2',
                    'tcella': 'tcell',
                    'tcellb': 'tcell2',
                    }
compartments = good_compartments.union(list(bad_compartments.keys()))

IPI_REGEX = re.compile(r'^(?P<ipi>(IPI[A-Z]{2,4}[0-9]{3}))[.]'
                       r'(?P<tn>(([TNRML]|NASH)[0-9]))[.]'
                       r'rna[.]'
                       r'(?P<stain>(' + '|'.join(compartments) + '))'
                       r'(?P<plural>(s*))'
                       r'(?P<puro>([Pp]*))'
                       r'(?P<suffix>([1-9]*))$')

data = {}
for tpm_file in glob.glob(os.path.join(os.environ['KLAB'], 'arrao/data/gene_expression/tpms/*')):
    plate = os.path.basename(tpm_file).split('_')[0]
    data[plate] = list(pd.read_csv(tpm_file, header=0, index_col=0, sep='\t', nrows=0).columns)
    assert len(set(data[plate])) == len(data[plate])


fix_commands = []
for plate in data:
    for sample in data[plate]:
        sample_re = (IPI_REGEX.search(sample) or 
                  CTRL_REGEX.search(sample) or 
                  BAD_CTRL_REGEX.search(sample))
        if not sample_re:
            if sample == 'UHR_CONTROL.plate4':
                fix_commands.append(rename_script('plate4', 'UHR_CONTROL.plate4', 'CONTROL_UHR.plate4'))
                continue
            raise RuntimeError
        if not sample_re.groupdict():
            # good control regex
            continue
        else:
            fail = False
            if 'ctrl' in sample_re.groupdict():
                fail = True
                fixed_name = 'CONTROL_{}.Plate{}'.format(
                    'UHR' if sample_re['juhr'].startswith(('U', 'u')) else 'jurkat',
                    sample_re['pltnum'] if sample_re['pltnum'].isnumeric() else sample_re['pltnum'].capitalize()
                    )
            else:
                stain = sample_re['stain']
                if sample_re['stain'] in bad_compartments:
                    stain = bad_compartments[sample_re['stain']]
                    fail = True
                if sample_re['plural']:
                    # Don't do anything.... the stain already is missing the's'
                    fail = True

                if sample_re['puro'] and sample_re['puro'].isupper():
                    fail=True
                    stain += 'p'
                
                if sample_re['suffix']:
                    # Add the suffix if it's greater than 1.
                    assert int(sample_re['suffix']) > 0
                    if int(sample_re['suffix']) > 1:
                        stain += sample_re['suffix']
                    # If it equals one, then dropping it is the expected output
                    else:
                        fail=True
                if fail:
                    fixed_name = '{}.{}.rna.{}'.format(sample_re['ipi'], sample_re['tn'], stain)
            if fail:
                fix_commands.append(rename_script(plate, sample, fixed_name))
                continue

print(fix_commands)