import fcsparser
import os
import sys

fcs_file = sys.argv[1]
print('procesing {}'.format(fcs_file))
fcs_file = os.path.abspath(fcs_file)
if not os.path.exists(fcs_file):
    raise RuntimeError('{} does not exist'.format(fcs_file))

try:
    meta, data = fcsparser.parse(fcs_file, reformat_meta=True)
except:
    print('Could not read {}'.format(fcs_file))
    raise()
else:
    for k in meta:
        print(k, meta[k], sep=': ')