import argparse
import glob
import gzip
import os
import pandas as pd
import pysam


def process_cel(input_prefix, output_prefix, chroms):
    print('Merging CEL files')
    final_cel_table = None
    barcodes_by_chrom = pd.DataFrame()
    index_name = None
    for chrom in chroms:
        print('Processing chrom {}'.format(chrom))
        df = pd.read_csv('{}{}.cel.gz'.format(input_prefix, chrom),
                          sep='\t',
                          header=0,
                          index_col=0)

        if final_cel_table is None:
            final_cel_table = df
            index_name = df.index.name
            col_order = df.columns
        else:
            final_cel_table = final_cel_table.append(df, ignore_index=True)

        df = df[['BARCODE']].copy()
        df[chrom] = df.index
        df.index = df['BARCODE']
        df.drop('BARCODE', axis=1, inplace=True)
        barcodes_by_chrom = pd.merge(barcodes_by_chrom, df,
                                     how='outer',
                                     left_index=True,
                                     right_index=True)

    barcodes_by_chrom = barcodes_by_chrom.fillna(-1).astype(int)
    barcodes_by_chrom['DROPLET_ID'] = range(len(barcodes_by_chrom.index))
    final_cel_table = final_cel_table.groupby('BARCODE').aggregate('sum')
    assert (final_cel_table.index == barcodes_by_chrom.index).all()
    final_cel_table['BARCODE'] = final_cel_table.index
    final_cel_table.index = barcodes_by_chrom['DROPLET_ID']
    final_cel_table.index.name = index_name
    final_cel_table = final_cel_table[col_order]

    final_cel_table.to_csv('{}.cel.gz'.format(output_prefix),
                           sep='\t',
                           header=True,
                           index=True)

    barcode_correction = {}
    for chrom in chroms:
        temp = barcodes_by_chrom.loc[barcodes_by_chrom[chrom]!=-1,[chrom, 'DROPLET_ID']].sort_values(by=chrom)
        barcode_correction[chrom] = dict(zip(temp[chrom], temp['DROPLET_ID']))

    return barcode_correction


def process_umi(input_prefix, output_prefix, chroms, barcode_correction):
    print('Merging UMI files')
    final_umi_table = pd.DataFrame(columns=['col1', 'col2', 'others'])
    i = 0
    for chrom in chroms:
        print('Processing chrom {}'.format(chrom))
        with gzip.open('{}{}.umi.gz'.format(input_prefix, chrom), 'rt') as iff:
            for line in iff:
                # Specifically nor strip()-ping
                line = line.split('\t')
                final_umi_table.loc[i] = [barcode_correction[chrom][int(line[0])], 
                                          line[1], 
                                          '\t'.join(line[2:])]
                i += 1

    final_umi_table.sort_values(['col1', 'col2'], inplace=True)
    with gzip.open('{}.umi.gz'.format(output_prefix), 'wt') as off:
        for x in  final_umi_table.itertuples():
            off.write('\t'.join([str(x.col1), x.col2, x.others]))

    return None


def process_var(input_prefix, output_prefix, chroms, barcode_correction):
    print('Merging VAR files')
    final_var_table = pd.DataFrame()
    i = 0
    var_correction = {}
    for chrom in chroms:
        print('Processing chrom {}'.format(chrom))
        df = pd.read_csv('{}{}.var.gz'.format(input_prefix, chrom),
                         sep='\t',
                         header=0,
                         index_col=0)
        
        var_correction[chrom] = dict(zip(df.index, 
                                         range(final_var_table.shape[0], 
                                               final_var_table.shape[0] + df.shape[0])))
        
        if final_var_table.empty:
            index_name = df.index.name
        final_var_table = final_var_table.append(df, ignore_index=True)


    final_var_table.index.name = index_name
    final_var_table.to_csv('{}.var.gz'.format(output_prefix),
                           sep='\t',
                           header=True,
                           index=True,
                           float_format="%.5f")

    return var_correction


def process_plp(input_prefix, output_prefix, chroms, barcode_correction, var_correction):
    print('Merging PLP files')
    final_plp_table = pd.DataFrame()
    
    for chrom in chroms:
        print('Processing chrom {}'.format(chrom))
        df = pd.read_csv('{}{}.plp.gz'.format(input_prefix, chrom),
                         sep='\t',
                         header=0,
                         index_col=0)
        df.index = df.index.map(lambda x: barcode_correction[chrom][x])
        df['SNP_ID'] = df['SNP_ID'].apply(lambda x: var_correction[chrom][x])
        if final_plp_table.empty:
            index_name = df.index.name
        final_plp_table = final_plp_table.append(df, ignore_index=False)

    final_plp_table.index.name = index_name
    final_plp_table['idx'] = final_plp_table.index
    final_plp_table.sort_values(['SNP_ID', 'idx'], inplace=True)
    final_plp_table.drop('idx', axis=1, inplace=True)

    final_plp_table.to_csv('{}.plp.gz'.format(output_prefix),
                           sep='\t',
                           header=True,
                           index=True)

    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_prefix')
    parser.add_argument('output_prefix')
    parser.add_argument('possorted_genome_bam')
    params = parser.parse_args()
    #params = parser.parse_args(['/scratch/arrao/foo/possorted_genome_bam_chr',
    #                            '/scratch/arrao/foo/possorted_genome_bam_MERGED',
    #                            '/krummellab/data1/immunox/MVIR1/data/single_cell_GEX/processed/MVIR1-POOL-SCG1/possorted_genome_bam.bam'])

    samfile = pysam.AlignmentFile(params.possorted_genome_bam, "rb")
    try:
        sorted_chroms = samfile.header.references
    finally:
        samfile.close()

    chroms = [c for c in sorted_chroms
                  if os.path.exists('{}{}.cel.gz'.format(params.input_prefix, c))]

    barcode_correction = process_cel(params.input_prefix, params.output_prefix, chroms)
    process_umi(params.input_prefix, params.output_prefix, chroms, barcode_correction)
    var_correction = process_var(params.input_prefix, params.output_prefix, chroms, barcode_correction)
    process_plp(params.input_prefix, params.output_prefix, chroms, barcode_correction, var_correction)
    

if __name__ == '__main__':
    main()

