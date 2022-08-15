import argparse
import gzip
import os
import pandas as pd
import shutil
import difflib

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freemuxlet_dir_tsv', required=True, help='A tsv file containing `sample`, '
                        '`freemuxlet_dir`, and optionally `barcode_suffix`. The suffix need not '
                        'contain field separator (e.g. _2) since we will automatically apply `--`.')
    parser.add_argument('--outdir', required=True, help='An output directory')
    params = parser.parse_args()

    fmx_ids = pd.read_table(params.freemuxlet_dir_tsv, sep="\t", index_col=None, header=0)
    if not {'sample', 'freemuxlet_dir'}.issubset(set(fmx_ids.columns)):
        raise RuntimeError('freemuxlet_dir_tsv must necessarily have columns named `sample` and '
                           '`freemuxlet_dir`')
    if 'barcode_suffix' not in fmx_ids.columns:
        if any(['--' in s for s in fmx_ids['sample']]):
            raise RuntimeError('samples inherently contain the `--` characters. Please explicitly '
                               'specify `barcode_suffix` in your tsv file.')
        fmx_ids['barcode_suffix'] = fmx_ids['sample']
    else:
        if any([s.startswith('-') for s in fmx_ids['barcode_suffix']]):
            raise RuntimeError('Values in `barcode_suffix` cannot begin with a `-`. We recommend '
                               'not using a field-separator at all.')

    if len(fmx_ids['barcode_suffix']) != len(set(fmx_ids['barcode_suffix'])):
            raise RuntimeError('Cannot have duplicate barcode suffixes. If freemuxlet_dir_tsv did '
                               ' not contain a column named `barcode_suffix`, you have duplicated '
                               '`sample` values and need to explicitly provide `barcode_suffix`.')
    
    if os.path.exists(params.outdir):
        outdir = os.path.join(params.outdir, 'merged_freemuxlet')
    else:
        assert os.path.exists(os.path.dirname(params.outdir))
        outdir = params.outdir
    os.mkdir(outdir)

    print('Processing *.var.gz')
    _x = None
    num_lines = {}
    # Collect nrows
    for s in fmx_ids.itertuples():
        with gzip.open(os.path.join(s.freemuxlet_dir, f'{s.sample}.var.gz')) as iff:
            x = len(iff.readlines())
        num_lines[s.sample] = x
    # If not all the same, check if smaller ones = largest just minus some lines at the end.  If so, that's okay!  Seems to be what dscpileup does when those SNPs are never captured.
    if len(set(list(num_lines.values()))) != 1:
        print(f'Found var.gz files of different lengths.')
        for x, y in num_lines.items():
            print(f'{x}: {y}')
        longest = max(num_lines, key=num_lines.get)
        with gzip.open(os.path.join(
            fmx_ids.loc[fmx_ids['sample']==longest, 'freemuxlet_dir'][0],
            f'{longest}.var.gz')
        ) as longest_file:
            longest_file_text = longest_file.readlines()
        for s in fmx_ids.itertuples():
            with gzip.open(os.path.join(s.freemuxlet_dir, f'{s.sample}.var.gz')) as this_file:
                this_file_text = this_file.readlines()
            diff = difflib.unified_diff(
                this_file_text,
                longest_file_text,
                fromfile=f'{s.sample}.var.gz',
                tofile=f'{longest}.var.gz'
            )
            print(diff)
            if len(diff) != 0 and len(diff) != num_lines[longest]-num_lines[s.sample]:
                print(diff)
                raise RuntimeError('Freemuxlet runs used different VCFs for processing')
    print('Writing merged.var.gz to disk')
    shutil.copy(os.path.join(s.freemuxlet_dir, f'{s.sample}.var.gz'),
                os.path.join(outdir, 'merged.var.gz'))

    merged_cels = None
    merged_plps = None
    print('Processing *.cel.gz and *.plp.gz concurrently')

    existing_cells = 0
    for s in fmx_ids.itertuples():
        print(f'Processing {s.sample}')
        cel = pd.read_csv(os.path.join(s.freemuxlet_dir, f'{s.sample}.cel.gz'), sep="\t", header=0,
                          index_col=None)
        plp = pd.read_csv(os.path.join(s.freemuxlet_dir, f'{s.sample}.plp.gz'), sep="\t", header=0,
                          index_col=None, dtype={
                              '#DROPLET_ID': int,
                              'SNP_ID': int,
                              'ALLELES': str,
                              'BASEQS': str,
                              })
        cel['#DROPLET_ID'] = cel['#DROPLET_ID'] + existing_cells
        plp['#DROPLET_ID'] = plp['#DROPLET_ID'] + existing_cells
        existing_cells += cel.shape[0]
        cel['BARCODE'] = cel['BARCODE'].apply(lambda x: f'{x}--{s.barcode_suffix}')

        if merged_cels is None:
            merged_cels = cel
            merged_plps = plp
        else:
            merged_cels = pd.concat([merged_cels, cel])
            merged_plps = pd.concat([merged_plps, plp])

    print('Writing merged.cel.gz to disk')
    merged_cels.to_csv(os.path.join(outdir, 'merged.cel.gz'),
                       header=True,
                       index=False,
                       sep="\t")
    print('Writing merged.barcodes.gz to disk')
    merged_cels.to_csv(os.path.join(outdir, 'merged.barcodes.gz'),
                       columns=['BARCODE'],
                       header=False,
                       index=False,
                       sep="\t")
    
    print('Sorting merged.plp.gz')
    merged_plps.sort_values(axis=0, by=['SNP_ID', '#DROPLET_ID'], inplace=True)
    print('writing merged.plp.gz to disk')
    merged_plps.to_csv(os.path.join(outdir, 'merged.plp.gz'), 
                       header=True,
                       index=False,
                       sep="\t")
    print('Done.')

if __name__ == '__main__':
    main()
