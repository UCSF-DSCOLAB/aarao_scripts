#!/bin/env python3
import argparse
import os
import screed

from time import time
from collections import Counter

def main():
    start_time = time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fq', type=str, help='Path to input Read 1', required=True)
    parser.add_argument('--error_threshold', type=float, help='Minimum value for fraction of most '
                        ' common barcode.', required=False, default=0.9)
    parser.add_argument('--num_read_to_process', type=int, help='Number of reads to process', 
                        required=False, default=1_000_000)
    params = parser.parse_args()

    assert os.path.exists(params.input_fq)
    assert params.num_read_to_process > 100_000, "Need at least 100k reads for good distribution"
    assert 0 < params.error_threshold <= 1 , "`error_threshold` must be in the range (0, 1] "

    barcode_counter = Counter()
    num_processed_reads = 0
    fq_file = screed.open(params.input_fq)
    try:
        for r in fq_file.iter_fn:
            barcode_counter[r['name'].split()[1].split(':')[-1]] += 1
            num_processed_reads += 1
            if num_processed_reads == params.num_read_to_process:
                break
    finally:
        fq_file.close()

    print(f'Processed {num_processed_reads} reads. Saw a total of {len(barcode_counter)} barcodes.',
          'Top 5 barcodes seen are:',
          sep='\n')

    for k, v in barcode_counter.most_common(5):
        print(f'\t{k}\t{v} ({round(100 * v/num_processed_reads,2)})%')

    if barcode_counter.most_common(1)[0][1] < params.error_threshold * num_processed_reads:
        print(f'WARNING: Most common barcode was < {100 * params.error_threshold}% of processed reads.')
          
    time_taken = time() - start_time
    print(f'Process took {int(time_taken//3600)}H:{int((time_taken//60)%60)}M:{int(time_taken%60)}S.',
          sep='\n')


if __name__ == '__main__':
    main()