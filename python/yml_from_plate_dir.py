import argparse
import os
import re
import glob

from collections import defaultdict

CONTROL_REGEX = re.compile(r'(?P<juhr>(control))', re.IGNORECASE)
JUHR_REGEX = re.compile(r'(?P<juhr>(jurkat|uhr))', re.IGNORECASE)
compartments = {'live',
                'tcell',
                'treg',
                'myeloid',
                'tumor',
                'stroma',
                'epcam',
                'microglia',
                'total',
                'cd11bposhladrneg',
                'cd45pos',
                'cd45',
                'dpost',
                'tpost',
                'tpre',
                'cd45neg',
                'teff',
                'tcellp',
                'livep'}

#IPI_REGEX = re.compile(r'(?P<ipi>(IPI[A-Z]{2,4}[0-9]{3}))')
#TNRML_REGEX = re.compile(r'(?P<tnrml>(([TNRML]|NASH)[0-9]))')
#COMPARTMENT_REGEX = re.compile(r'(?P<compartment>(' + '|'.join(compartments) + ')[s]*[2-9]*)')

IPI_REGEX = re.compile(r'(?P<ipi>(IPI[A-Z]{2,4}[0-9]{3}))_'
                       r'(?P<tnrml>(([TNRML]|NASH)[0-9]))_'
                       r'(?P<compartment>(' + '|'.join(compartments) + ')s*[2-9]*)_.*fastq.gz')


def parse_dir(input_dir, plate_num):
    all_files = glob.glob('{}/*.fastq.gz'.format(input_dir)) + glob.glob('{}/*/*.fastq.gz'.format(input_dir))
    results = defaultdict(list)
    for f in all_files:
        if CONTROL_REGEX.search(f):
            juhr = JUHR_REGEX.search(f)
            assert juhr is not None
            # Control sample
            results['CONTROL_{}.Plate{}'.format(juhr['juhr'].lower(), plate_num)].append(f)
        else:
            ipi = IPI_REGEX.search(f)
            assert ipi is not None 
            results['{}.{}.rna.{}'.format(ipi['ipi'].upper(), ipi['tnrml'].upper(), ipi['compartment'].lower())].append(f)
    assert {len(x) for x in results.values()}  == {2}
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('plate_num', help='plate num')
    parser.add_argument('input_dir', help='Input dir')
    parser.add_argument('--output_yml', help='Path to an output file', default=None)
    parser.add_argument('--overwrite', help='Overwrite existing yml?', action='store_true')
    params = parser.parse_args()

    assert params.plate_num in [str(x+1) for x in range(29)] + ['Nova']
    outfile = params.output_yml
    if outfile is not None:
        outfile = os.path.abspath(params.output_yml)
        if os.path.exists(outfile) and not params.overwrite:
            raise RuntimeError('Cowardly refusing to overwrite {}'.format(outfile))

    input_dir = os.path.abspath(params.input_dir)
    assert os.path.exists(params.input_dir)

    files = parse_dir(input_dir, params.plate_num)
    off = open(outfile, 'w') if outfile else None
    try:
        for sample, fastqs in files.items():
            print('- :sample_name: {}'.format(sample), file=off)
            print('  :sample_files:', file=off)
            for f in sorted(fastqs):
                print('    - {}'.format(f), file=off)
    except:
        if outfile:
            outfile.close()
                
if __name__ == '__main__':
    main()
