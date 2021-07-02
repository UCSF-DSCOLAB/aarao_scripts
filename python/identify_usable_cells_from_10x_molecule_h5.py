import argparse
import gzip
import os
import tables
from collections import defaultdict, Counter

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cellranger_molecule_h5', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--GEX_umi_cutoff', type=int, default=0, required=False)
    parser.add_argument('--GEX_feature_cutoff', type=int, default=100, required=False)
    parser.add_argument('--overwrite', action-'store_true')
    params = parser.parse_args()

    cellranger_molecule_h5 = os.path.abspath(params.cellranger_molecule_h5)
    assert os.path.exists(cellranger_molecule_h5) and os.path.isfile(cellranger_molecule_h5)

    outdir = os.path.abspath(outdir)
    assert os.path.exists(outdir) and os.path.isdir(outdir)

    assert (params.GEX_umi_cutoff >=0 and 
            params.GEX_feature_cutoff >=0 and 
            params.GEX_umi_cutoff + params.GEX_feature_cutoff > 0)

    outfile = os.path.join(outdir, 'barcode_of_interest.list.gz')
    if os.path.exists(outfile):
        if params.overwrite:
            raise RuntimeError('ERROR: Outfile {} exists in outdir {}. Cannot overwrite unless '
                               '--overwrite is specified'.format(os.path.basename(outfile), outdir))
        else:
            print('WARNING: overwriting existing outfile {}'.format(outfile))

    with tables.open_file(cellranger_molecule_h5, 'r') as iff:
        try:
            barcodes = iff.get_node(iff.root, 'barcodes')
            barcode_idx = iff.get_node(iff.root, 'barcode_idx')
            umi_counts = iff.get_node(iff.root, 'count')
            feature_idx = iff.get_node(iff.root, 'feature_idx')
            feature_type = iff.get_node(iff.root, 'features/feature_type')

        except tables.NoSuchNodeError:
            print "Matrix group does not exist in this file."
            return None

        per_bc_umi_counts = Counter()
        per_bc_feature_counts = defaultdict(set)
        for i in range(0, barcode_idx.shape[0], 1):
            if feature_type[feature_idx[i]] != b'Gene Expression':
                continue            
            per_bc_umi_counts[barcode_idx[i]] += umi_counts[i]
            per_bc_feature_counts[barcode_idx[i]].add(feature_idx[i])

        per_bc_feature_counts = {b: len(v) for b, v in per_bc_feature_counts.items()}
        
        barcode_idxs_of_interest_by_umi = {b for b, v in per_bc_umi_counts.items() 
                                             if v >= params.GEX_umi_cutoff}
        barcode_idxs_of_interest_by_feature = {b for b, v in per_bc_feature_counts.items() 
                                                 if v >= params.GEX_feature_cutoff}
        
        barcode_idxs_of_interest = 
            list(barcode_idxs_of_interest_by_umi & barcode_idxs_of_interest_by_feature)

        barcodes_of_interest = barcodes[barcode_idxs_of_interest]
        with gzip.open(outfile, 'wb') as off:
            for barcode in barcodes_of_interest:
                off.write(barcode + b'\n')

if __name__ = ='__main__':
    main()