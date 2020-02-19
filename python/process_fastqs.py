from __future__ import print_function

import argparse
import glob
import logging
import multiprocessing as mp
import os
import subprocess
import sys

from functools import partial

def process_fastq(worker_id, fastqs, results):
    print('INFO: FASTQ processing worker %s is up and running.' % worker_id)
    while True:
        fastq_file = fastqs.get(timeout=30)
        if fastq_file is None:
            break
        
        try:
            md5sum = subprocess.check_output(['md5sum',fastq_file])
        except subprocess.CalledProcessError:
            results[os.path.basename(fastq_file)] = (None, False)
            continue

        md5sum = md5sum.split()[0]

        try:
            integrity_ok = subprocess.check_output(['gzip', '--test',fastq_file])
        except subprocess.CalledProcessError:
            results[os.path.basename(fastq_file)] = (md5sum, False)
            continue
        
        results[os.path.basename(fastq_file)] = (md5sum, True)

    print('INFO: FASTQ processing worker %s received signal to go down.' % worker_id)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('plate', help='plate number to process')
    parser.add_argument('--fastq_dir', help='Directory containing all fastqs.',
                        default='/cbc2/data2/samadb/IPI/fastqs', required=False)
    parser.add_argument('--output_dir', help='Directory containing output md5 sums.',
                        default='/data/shared/krummellab/ipi/md5sums', required=False)
    
    params = parser.parse_args()

    manager = mp.Manager()
    fastqs = manager.Queue()
    results = manager.dict()

    fastq_path = os.path.join(params.fastq_dir, 'Plate%sRnaseq' % params.plate)

    assert os.path.exists(fastq_path), '%s does not point to a valid Directory' % fastq_path

    if 'PBS_NUM_PPN' in os.environ:
        cores = int(os.environ['PBS_NUM_PPN'])
    else:
        cores = mp.cpu_count()
    print('INFO: Identified %s cpu cores.' % cores)

    for fastq_file in glob.glob(os.path.join(fastq_path,'*','*.gz')):
        fastqs.put(fastq_file)

    for _ in range(0, cores):
        # Add the graceful worker shutdown signal to the queue for every worker.
        fastqs.put(None)

    print('INFO: Starting %s workers for fastq processing.' % cores)
    
    pool = mp.Pool(processes=cores)
    process_fastq_partial = partial(process_fastq, fastqs=fastqs, results=results)
    pool.map(process_fastq_partial, range(0, cores))
    pool.close()
    pool.join()

    with open(os.path.join(params.output_dir, 'Plate_%s_md5sums.txt' % params.plate), 'w') as off:
        for fastq_file in results.keys():
            md5, integrity_ok = results[fastq_file]
            print(fastq_file, md5 if integrity_ok else 'ERROR', sep='\t', file=off)

    print('INFO: Completed. Processed %s fastq files.' % len(results))

if __name__ == '__main__':
    main()
