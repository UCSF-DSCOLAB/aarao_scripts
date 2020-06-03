import argparse
import gzip
import os
import pandas as pd
import pysam

from collections import Counter


def count_umis(bamfile, bedfile, interval_key, outfile, read_update_verbosity=100000, 
               interval_update_verbosity=100, mapq_threshold=255, min_overlap_bases=0):
    
    intervals = pd.read_csv(bedfile, sep='\t', 
                            header=None, 
                            names=['chrom', 
                                   'start',
                                   'end'], 
                            dtype={'chrom': str, 
                                   'start': int,
                                   'end': int})

    udict = {}
    sam_handle = pysam.Samfile(bamfile, 'rb')

    noCB_readnum = 0
    bad_umi = 0
    umi_dup = 0
    good_read = 0
    num_cells = 0
    duplicate = 0
    map_qual = Counter()
    low_overlap = 0

    reads_processed = 0
    intervals_processed = 0
    try:
        for interval in intervals.itertuples():
            intervals_processed += 1
            if intervals_processed % (max(50, interval_update_verbosity/100)) == 0:
                print('Processed {} intervals'.format(intervals_processed))
            if (interval.end - interval.start) <  min_overlap_bases:
                continue
                
            for read in sam_handle.fetch(interval.chrom, start=interval.start, stop=interval.end):
                reads_processed += 1
                
                if reads_processed % read_update_verbosity == 0:
                    print('Processed {} reads'.format(reads_processed))
                if read.is_duplicate:
                    duplicate += 1
                    continue

                map_qual[read.mapping_quality] += 1
                
                if read.mapping_quality < mapq_threshold:
                    continue

                if sum([interval.start <= r <=interval.end 
                            for r in read.get_reference_positions()]) < min_overlap_bases:
                    low_overlap += 1
                    continue 

                rdict = dict(read.tags)
                if 'CB' not in rdict:
                    noCB_readnum += 1
                else:
                    if 'UB' not in rdict:
                        # Bad Umi... skip
                        bad_umi += 1
                        continue
                    if rdict['CB'] in udict:
                        if rdict['UB'] in udict[rdict['CB']]:
                            # This is a umi duplicate
                            umi_dup += 1
                            continue
                        else:
                            good_read += 1
                            udict[rdict['CB']].add(rdict['UB'])
                    else:
                        good_read += 1
                        num_cells += 1
                        udict[rdict['CB']] = {rdict['UB']}
    except:
        print('Failed on read {}'.format(read.query_name))
        print(read.to_dict())
        print('######')
        print(read.tags)
        raise
    finally:
        sam_handle.close()
        print('Identified:',              
              'Unique Cells:{}'.format(num_cells),
              
              'Total reads processed: {}'.format(reads_processed),
              '    Good reads: {}'.format(good_read),
              '    Duplicated reads: {}'.format(duplicate),
              '    Reads rejected for low overlap with region: {}'.format(low_overlap),
              '    Reads arising from a cell without a proper barcodes: {}'.format(noCB_readnum),
              '    Reads with bad UMIs: {}'.format(bad_umi),
              '    UMI duplicates: {}'.format(umi_dup), 
              '    Reads with low map qual: {}'.format(sum([v for m, v in map_qual.items() if m < mapq_threshold])),
              '    MAPQ_distribution:\n        {}'.format('\n        '.join(['{} : {}'.format(m, v) for m, v in map_qual.items()])),
              sep='\n')
    cdict = {x: len(y) for x,y in udict.items()}
    
    with open(outfile, 'w') as off:
        print('barcode', interval_key, sep='\t', file=off)
        for bc, cnt in cdict.items():
            print(bc, cnt, sep='\t', file=off)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_bam', required=True)
    parser.add_argument('--intervals_bed', help="A  bed file with 3 columns per line, chrom, start "
                        ", end. Must be 0-based indexing.", required=True)
    parser.add_argument('--interval_key', help="The key to use for the output metadata tsv. "
                        "Inferred from the name `bedfile` if not provided", required=False, 
                        default="infer")
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--overwrite', action="store_true")
    parser.add_argument('--read_update_verbosity', type=int, default=1000000)
    parser.add_argument('--interval_update_verbosity', type=int, default=1000)
    parser.add_argument('--mapq_threshold', type=int, default=255)
    parser.add_argument('--min_overlap_bases', type=int, default=0)
    
    params = parser.parse_args()

    input_bam = os.path.abspath(params.input_bam)
    assert os.path.exists(input_bam)
    
    intervals_bed = os.path.abspath(params.intervals_bed)
    assert os.path.exists(intervals_bed)

    if params.interval_key == 'infer':
        interval_key = os.path.splitext(os.path.basename(intervals_bed))[0]
    else:
        interval_key = params.interval_key

    outfile = os.path.abspath(params.outfile)
    if os.path.exists(params.outfile):
        if params.overwrite:
            print('Overwriting existing file: {}'.format(outfile))
        else:
            raise RuntimeError('File exists ({}). specify --overwrite to overwrite '.format(outfile))

    count_umis(input_bam, intervals_bed, interval_key, outfile, params.read_update_verbosity,
               params.interval_update_verbosity, params.mapq_threshold, params.min_overlap_bases)

if __name__ == '__main__':
    main()