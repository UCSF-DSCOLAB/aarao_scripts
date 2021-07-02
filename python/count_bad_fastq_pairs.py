#!/bin/env python3
import argparse
import gzip
import os
import re
import screed
import tempfile

from time import time
from collections import OrderedDict

class PairedScreedFastqs(object):
    def __init__(self, seqfile1, seqfile2):
        assert os.path.exists(seqfile1)
        assert os.path.exists(seqfile2)
        self.seqfile1 = screed.open(seqfile1)
        self.seqfile2 = screed.open(seqfile2)
        self.balanced = True
    def __iter__(self): 
        return self
    def __next__(self):
        while True:
            exhausted = None
            try:
                r1 = next(self.seqfile1.iter_fn)
            except StopIteration:
                exhausted = 1
                r1 = None
            try:
                r2 = next(self.seqfile2.iter_fn)
            except StopIteration:
                exhausted = 2
                r2 = None
            if r1 is None and r2 is None:
                raise StopIteration
            elif exhausted is not None:
                print(f'WARNING: Reached EOF in R{exhausted}.')
                self.balanced = False
            return r1, r2
    def close(self):
        self.seqfile1.close()
        self.seqfile2.close()
    def __del__(self):
        self.close()

def count_bad_fastq_pairs(in_fq1, in_fq2, update_interval=1_000_000):
    # Prepare the input fastq generator
    fqpairs = PairedScreedFastqs(in_fq1,
                                 in_fq2)

    # Logging variables
    bad_pairs = _bad_pairs = processed_reads = 0
    last_update_time = time()

    # Wrap the rest in a try: finally struct so we can close files appropriately
    try:
        for r1, r2 in fqpairs:
            if processed_reads % update_interval == 0:
                print(f'\tProcessed {processed_reads} read pairs '
                      f'({round(60 * update_interval/(time()-last_update_time), 2)} pairs/min)... '
                      f'Identified {_bad_pairs} pairs since last update.')
                bad_pairs += _bad_pairs
                _bad_pairs = 0
            processed_reads += 1
            r1_name = r1['name'].split()[0] if r1 else None
            r2_name = r2['name'].split()[0] if r2 else None
            # Best case scenario. Properly paired!
            _bad_pairs += (r1_name != r2_name)
        bad_pairs += _bad_pairs
    finally:
        fqpairs.close()
    return processed_reads, bad_pairs, fqpairs.balanced


def main():
    start_time = time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fq1', type=str, help='Path to input Read 1', required=True)
    parser.add_argument('-I', '--input_fq2', type=str, help='Path to input Read 2', required=True)
    parser.add_argument('--update_interval', type=int, help='Number of read pairs processed per update', 
                        required=False, default=1_000_000)
    params = parser.parse_args()

    assert os.path.exists(params.input_fq1)
    assert os.path.exists(params.input_fq2)
    processed_reads, bad_pairs, balanced = count_bad_fastq_pairs(params.input_fq1, params.input_fq2, 
                                                                 update_interval=params.update_interval)
    
    print(f'Completed Process....',
          f'Processed {processed_reads} read pairs and identified {bad_pairs} bad pairs.',
          sep='\n')
    if bad_pairs > 0:
        print('ERROR: Identified bad pairing in the provided fastq set.')
    if not balanced:
        print('ERROR: The input read files are not balanced.')

    time_taken = time() - start_time
    print(f'Process took {int(time_taken//3600)}H:{int((time_taken//60)%60)}M:{int(time_taken%60)}S.')
    return bad_pairs != 0 or not balanced

if __name__ == '__main__':
    main()
