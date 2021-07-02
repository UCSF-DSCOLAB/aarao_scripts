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
            return r1, r2
    def close(self):
        self.seqfile1.close()
        self.seqfile2.close()
    def __del__(self):
        self.close()


def write_pair_to_files (r1, r2, o1, o2):
    if r1 is not None:
        write_record_to_file(r1, o1)
    if r2 is not None:
        write_record_to_file(r2, o2)


def write_record_to_file(r, o):
    print("@{name}\n{sequence}\n+{annotations}\n{quality}".format(**r), file=o)


# This is a recursive function
def parse_fastq_pair(in_fq1, in_fq2, out_fq1, out_fq2, sngltn_fq1, sngltn_fq2, 
                     max_reads_in_memory=200_000, update_interval=1_000_000, 
                     cache_ejection_rate=0.1, tempdir=None, delete_input_fastqs=False, iter_num=1):
    start_time = time()
    # Prepare the input fastq generator
    fqpairs = PairedScreedFastqs(in_fq1,
                                 in_fq2)
    # Open the descriptors for the properly paired output files
    off1 = gzip.open(out_fq1, 'at')
    off2 = gzip.open(out_fq2, 'at')
    # Open the descriptors for the Singleton files
    sff1 = gzip.open(sngltn_fq1, 'at')
    sff2 = gzip.open(sngltn_fq2, 'at')
    # Create temporary spillover files
    soff1 = tempfile.mkstemp(prefix='pair_fastqs', dir=tempdir)
    os.close(soff1[0])  # To prevent unintentional corruption from GC
    soff1 = gzip.open(soff1[1], 'at')
    soff2 = tempfile.mkstemp(prefix='pair_fastqs', dir=tempdir)
    os.close(soff2[0])  # To prevent unintentional corruption from GC
    soff2 = gzip.open(soff2[1], 'at')
    # Ordered Dicts of unpaired reads so we can spill to disk in a FIFO manner
    r1_reads_in_mem = OrderedDict()
    r2_reads_in_mem = OrderedDict()
    # Store names of reads spilled over to disk so we can identify singletons early
    spilled_over_r1_readnames = set()
    spilled_over_r2_readnames = set()
    # Logging variables
    written_pairs = processed_reads = 0
    spilled_over_r1_num = spilled_over_r2_num = 0
    singleton_r1_num = singleton_r2_num = 0
    last_update_time = time()
    # Wrap the rest in a try: finally struct so we can close files appropriately
    try:
        for r1, r2 in fqpairs:
            processed_reads += 1
            if processed_reads % update_interval == 0:
                print(f'\tProcessed {processed_reads} read pairs '
                      f'({round(60 * update_interval/(time()-last_update_time), 2)} pairs/min)... '
                      f'wrote {written_pairs} to disk and cache size is currently '
                      f'{len(r1_reads_in_mem) + len(r2_reads_in_mem)} reads')
            r1_name = r1['name'].split()[0] if r1 else None
            r2_name = r2['name'].split()[0] if r2 else None
            # Best case scenario. Properly paired!
            if r1_name == r2_name:
                written_pairs += 1
                write_pair_to_files(r1, r2, off1, off2)
                continue
            # If we got an R1 (FQ1 isn't exhausted)
            if r1_name is not None:
                # If the mate is in memory write it out else add R1 to memory
                if r1_name in r2_reads_in_mem:
                    written_pairs += 1
                    write_pair_to_files(r1, r2_reads_in_mem.pop(r1_name), off1, off2)
                else:
                    r1_reads_in_mem[r1_name] = r1
            # If we got an R2 (FQ2 isn't exhausted)
            if r2_name is not None:
                # If the mate is in memory write it out else add R2 to memory
                if r2_name in r1_reads_in_mem:
                    written_pairs += 1
                    write_pair_to_files(r1_reads_in_mem.pop(r2_name), r2, off1, off2)
                else:
                    r2_reads_in_mem[r2_name] = r2
            # Spillover logic
            if len(r1_reads_in_mem) + len(r2_reads_in_mem) >= max_reads_in_memory:
                print('Cache exceeded max size. Spilling over some recors to disk.')
                # Drop ~20% of all reads in the order they were added
                for _ in range(int(cache_ejection_rate * max_reads_in_memory)):
                    if len(r1_reads_in_mem) > 0:
                        k, v = r1_reads_in_mem.popitem(last=False)
                        spilled_over_r1_readnames.add(k)
                        write_record_to_file(v, soff1)
                    if len(r2_reads_in_mem) > 0:
                        k, v = r2_reads_in_mem.popitem(last=False)
                        spilled_over_r2_readnames.add(k)
                        write_record_to_file(v, soff2)
        # If we reach here, we've exhausted the input fastq files. Any reads in memory are not 
        # paired with values in the other read (else they would have been ejected), but may still be
        # paired with a read that was spilled over to disk. If they aren't, then write them out as a
        # Singleton now
        spilled_over_r1_num += len(spilled_over_r1_readnames)
        spilled_over_r2_num += len(spilled_over_r2_readnames)
        for k, v in r1_reads_in_mem.items():
            if k in spilled_over_r2_readnames:
                write_record_to_file(v, soff1)
                spilled_over_r1_num += 1
            else:
                write_record_to_file(v, sff1)
                singleton_r1_num += 1
        for k, v in r2_reads_in_mem.items():
            if k in spilled_over_r1_readnames:
                write_record_to_file(v, soff2)
                spilled_over_r2_num += 1
            else:
                write_record_to_file(v, sff2)
                singleton_r2_num += 1
        # All records are accounted for.
    finally:
        fqpairs.close()
        off1.close()
        off2.close()
        sff1.close()
        sff2.close()
        soff1.close()
        soff2.close()
        if delete_input_fastqs:
            os.remove(off1.name)
            os.remove(off2.name)
        time_taken = time() - start_time
        print(f'Completed Iteration {iter_num}...',
              f'Processed {processed_reads} read pairs.',
              f'\tWrote {written_pairs} pairs to the output file.',
              f'\tSpilled over {spilled_over_r1_num} reads from R1.',
              f'\tSpilled over {spilled_over_r2_num} reads from R2.',
              f'\tIdentified {singleton_r1_num} singletons from R1.',
              f'\tIdentified {singleton_r2_num} singletons from R2.',
              f'Process took {int(time_taken//3600)}H:{int((time_taken//60)%60)}M:{int(time_taken%60)}S.',
              sep='\n')
    if spilled_over_r1_num + spilled_over_r2_num > 0:
        if iter_num > 5:
            raise RuntimeError('Took 5 iters and still not done.... maybe rerun with a higher '
                               'value for `max_reads_in_memory`')
        # There's still reads left to process
        wp, sr1n, sr2n = parse_fastq_pair(soff1.name, soff2.name, 
                                          out_fq1, out_fq2, 
                                          sngltn_fq1, sngltn_fq2, 
                                          max_reads_in_memory=max_reads_in_memory, 
                                          update_interval=update_interval, tempdir=tempdir,
                                          delete_input_fastqs=True, iter_num=iter_num+1)
        written_pairs += wp
        singleton_r1_num += sr1n
        singleton_r2_num += sr2n
    else:
        # Base case
        # These files should be empty but clean them up anyway
        # The files won't be empty (they're empty but gzippped so they still are > 0 bytes)
        # However they should be the same size on a given machine
        assert os.stat(soff1.name).st_size == os.stat(soff2.name).st_size
        os.remove(soff1.name)
        os.remove(soff2.name)
    return written_pairs, singleton_r1_num, singleton_r2_num


def main():
    start_time = time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fq1', type=str, help='Path to input Read 1', required=True)
    parser.add_argument('-I', '--input_fq2', type=str, help='Path to input Read 2', required=True)
    parser.add_argument('--out_prefix', type=str, help='Path to a prefix for output and singleton '
                        'files. Will be overruled by values to -o/O/s/S', 
                        required=False, default='./properly_paired')
    parser.add_argument('-o', '--output_fq1', type=str, help='Path to output Read 1', 
                        required=False, default=None)
    parser.add_argument('-O', '--output_fq2', type=str, help='Path to output Read 2', 
                        required=False, default=None)
    parser.add_argument('-s', '--singleton_fq1', type=str, help='Path to singletons from Read 1', 
                        required=False, default=None)
    parser.add_argument('-S', '--singleton_fq2', type=str, help='Path to singletons from Read 2', 
                        required=False, default=None)
    parser.add_argument('--overwrite', action='store_true', help='Can we overwrite existing outputs?')
    parser.add_argument('--max_reads_in_memory', type=int, help='Maximum read to store in memory '
                        '(R1+R2) before spilling to disk', required=False, default=200_000)
    parser.add_argument('--cache_ejection_rate', type=float, help='Fraction of cache to spill to '
                        'disk if size gets too large.', required=False, default=0.1)
    parser.add_argument('--tempdir', type=str, help='Path to a tempdir to write excess reads', 
                        required=False, default=None)
    parser.add_argument('--update_interval', type=int, help='Number of read pairs processed per update', 
                        required=False, default=1_000_000)
    params = parser.parse_args()
    #params = parser.parse_args(['-i', 'test1.fq.gz', '-I', 'test2.fq.gz'])

    assert os.path.exists(params.input_fq1)
    assert os.path.exists(params.input_fq2)
    assert params.max_reads_in_memory >= 1000
    assert params.update_interval >= 1000
    assert 0 < params.cache_ejection_rate <= 1

    prefix = os.path.abspath(params.out_prefix)
    
    output_fq1 = os.path.abspath(params.output_fq1) if params.output_fq1 else f'{prefix}_out_1.fq.gz'
    output_fq2 = os.path.abspath(params.output_fq2) if params.output_fq2 else f'{prefix}_out_2.fq.gz'
    singleton_fq1 = os.path.abspath(params.singleton_fq1) if params.singleton_fq1 else f'{prefix}_singletons_1.fq.gz'
    singleton_fq2 = os.path.abspath(params.singleton_fq2) if params.singleton_fq2 else f'{prefix}_singletons_2.fq.gz'

    if any([os.path.exists(output_fq1), os.path.exists(output_fq2), 
            os.path.exists(singleton_fq1), os.path.exists(singleton_fq2)]):
        print('Some output paths exist....',
              f'\toutput_fq1 ({output_fq1}) exists: {os.path.exists(output_fq1)}', 
              f'\toutput_fq2 ({output_fq2}) exists: {os.path.exists(output_fq2)}', 
              f'\tsingleton_fq1 ({singleton_fq1}) exists: {os.path.exists(singleton_fq1)}', 
              f'\tsingleton_fq2 ({singleton_fq2}) exists: {os.path.exists(singleton_fq2)}',
              sep='\n')
        if not params.overwrite:
            raise RuntimeError('ERROR: `--overwrite` was not provided. Cannot continue.')
        else:
            print('WARNING: Deleting existing output files so they may be regenerated.')
            if os.path.exists(output_fq1):
                os.remove(output_fq1)
            if os.path.exists(output_fq2):
                os.remove(output_fq2)
            if os.path.exists(singleton_fq1):
                os.remove(singleton_fq1)
            if os.path.exists(singleton_fq2):
                os.remove(singleton_fq2)

    wp, sr1n,sr2n = parse_fastq_pair(params.input_fq1, params.input_fq2, 
                                     output_fq1, output_fq2, 
                                     singleton_fq1, singleton_fq2, 
                                     max_reads_in_memory=params.max_reads_in_memory, 
                                     update_interval=params.update_interval,
                                     cache_ejection_rate=params.cache_ejection_rate,
                                     tempdir=params.tempdir, 
                                     delete_input_fastqs=False, 
                                     iter_num=1)
    time_taken = time() - start_time
    print(f'Completed Process....',
          f'Processed a total of {wp * 2 + sr1n + sr2n} reads.',
              f'\tWrote {wp} pairs to the output file.',
              f'\tIdentified {sr1n} singletons from R1.',
              f'\tIdentified {sr2n} singletons from R2.',
              f'Process took {int(time_taken//3600)}H:{int((time_taken//60)%60)}M:{int(time_taken%60)}S.',
              sep='\n')


if __name__ == '__main__':
    main()