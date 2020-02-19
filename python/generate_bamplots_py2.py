from __future__ import division, print_function

import argparse
import gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import pysam

from collections import Counter
from functools import partial
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Manager, Pool



def sam_bam_cram(seqfile):
    """
    Is the input file sam (''), bam ('b') or cram ('c')
    """
    with open(seqfile) as iff:
         x = iff.read(4)
    if x == 'CRAM':
        return 'c'
    elif x[:3] == '\x1f\x8b\x08':
        # Could be a regular gzip file. Ensure it's actually a bam
        with gzip.open(seqfile, 'r') as iff:
            x = iff.read(4)
        if x == 'BAM\x01':
            return 'b'
        else:
            assert False, "%s is (b)gzipped but is not a bam" % seqfile
    elif x[:3] == "@HD":
        return ''
    else:
        assert False, "%s is neither sam, bam, or cram"  % seqfile


def _get_mapping_data_multi_threaded(worker_id, seqfile, reference_file, queue, map_stats, lock):
    print('Worker %s is up and running.' % worker_id)
    samfile = pysam.Samfile(seqfile, 'r' + sam_bam_cram(seqfile), reference_filename=reference_file)
    try:
        while True:
            line = queue.get(timeout=30)
            if line is None:
                break
            _map_stats = {'mapqs': Counter(), 'insert_sizes': Counter()}
            for read in samfile.fetch(contig=line):
                _map_stats['mapqs'][-1 if read.is_unmapped else read.mapping_quality] += 1
                if read.is_read2:
                    continue
                if read.is_proper_pair:
                    _map_stats['insert_sizes'][abs(read.isize)] += 1

            lock.acquire()
            for _dict in 'mapqs', 'insert_sizes':
                if _dict in map_stats.keys():
                    _map_stats[_dict].update(map_stats[_dict])
                    map_stats[_dict] = _map_stats[_dict]
                else:
                    map_stats[_dict] = _map_stats[_dict]
            lock.release()
            
    finally:
        samfile.close()
    print('Worker %s received signal to go down.' % worker_id)

def _get_mapping_data_single_threaded(seqfile, reference_file, indexed):
    samfile = pysam.Samfile(seqfile, 'r' + sam_bam_cram(seqfile), reference_filename=reference_file)
    map_stats = {'mapqs': Counter(), 'insert_sizes': Counter()}
    try:
        for read in samfile.fetch(until_eof=not indexed):
            map_stats['mapqs'][-1 if read.is_unmapped else read.mapping_quality] += 1
            if read.is_read2 or read.is_unmapped:
                continue
            if read.is_proper_pair:
                map_stats['insert_sizes'][abs(read.isize)] += 1
    finally:
        samfile.close()
    return map_stats


def get_mapping_data(seqfile, num_threads, reference_file=None):
    assert num_threads >= 1

    samfile = pysam.Samfile(seqfile, 'r' + sam_bam_cram(seqfile), reference_filename=reference_file)
    try:
        samfile.check_index()
    except ValueError as e:
        if e.message == 'mapping information not recorded in index or index not available':
            indexed = False
        else:
            raise
    else:
        indexed = True
    finally:
        samfile.close()
   
    if indexed and num_threads > 1:
        manager = Manager()
        queue = manager.Queue()
        map_stats = manager.dict()
        lock = manager.Lock()
        
        samfile = pysam.Samfile(seqfile, 'r' + sam_bam_cram(seqfile),
                                reference_filename=reference_file)
        for contig in samfile.references:
            queue.put(contig)
        for _ in range(0, num_threads):
            # Add the graceful worker shutdown signal to the queue for every worker.
            queue.put(None)
        samfile.close()
        pool = Pool(processes=num_threads)
        get_mapping_data_partial = partial(_get_mapping_data_multi_threaded,
                                           seqfile=seqfile,
                                           reference_file=reference_file,
                                           queue=queue,
                                           map_stats=map_stats,
                                           lock=lock)
        pool.map(get_mapping_data_partial, range(0, num_threads))
        pool.close()
        pool.join()
    else:
        map_stats = _get_mapping_data_single_threaded(seqfile, reference_file, indexed)

    return Counter(map_stats['mapqs']), Counter(map_stats['insert_sizes'])



def generate_mapplots(seqfile, num_threads=1, reference_file=None, skip_mapq=False,
                      skip_insert_sizes=False):
    assert not (skip_mapq and skip_insert_sizes)
    assert num_threads >0
    out_prefix = os.path.splitext(seqfile)[0]

    mapqs, insert_sizes = get_mapping_data(seqfile, num_threads, reference_file)

    unmapped = mapqs[-1]
    qmax_pct = round(100 * mapqs[max(mapqs)]/sum([mapqs[x] for x in mapqs if x != -1]), 2)

    with PdfPages(out_prefix + '_bamqc.pdf') as pdf:
        # First plot the mapq
        fig, ax = plt.subplots(figsize=(20,20))
        x, h = zip(*sorted(mapqs.items(), key=lambda z: z[0]))
        x = [str(_x) for _x in x]
        bar = ax.bar(x, h)
        ax.set_xlabel('MAPQ')
        ax.set_ylabel('Count')
        ax.text(bar[-1].get_x() + bar[-1].get_width()/2.0, bar[-1].get_height(), '%s%%' % qmax_pct, 
                horizontalalignment='center', verticalalignment='bottom', fontsize=15)
        plt.xticks(rotation=45)
        plt.title('MapQ profile')
        pdf.savefig(dpi=600)
        #plt.savefig(out_prefix + '_mapq.pdf', dpi=600)

        # Second plot the insert sizes
        fig, ax = plt.subplots(figsize=(20,10))
        x, y = zip(*sorted(insert_sizes.items(), key=lambda z: z[0]))
        ax.plot(x, y)
        ax.set_xlabel('Insert Size')
        ax.set_ylabel('Count')

        # An insert size of 0 means there was an alignment to another chromosome, or an unmapped 
        # read.  We should not consider this when estimating an insert size.
        if x[0] == 0:
            x = list(x)
            x.pop(0)
            y = list(y)
            y.pop(0)
        total_events = sum(insert_sizes.keys())
        
        y_sum = sum(y)
        median = [x[i] for i, j in enumerate(y) if sum(y[:i+1]) >= y_sum/2][0]
        ax.axvline(x=median, linestyle='--', color='red', linewidth=0.5)
        ax.text(median, 0, 'median=%s' % median, horizontalalignment='right',
                verticalalignment='top')

        _total_events = total_events
        percentile = 0.99
        for i in sorted(insert_sizes.keys(), reverse=True):
            _total_events -= insert_sizes[i]
            if _total_events/total_events <= percentile:
                break

        ax.set_xlim(0,i)
        plt.title('Insert Size profile')
        pdf.savefig(dpi=600)
        #plt.savefig(out_prefix + '_insert_sizes.pdf', dpi=600)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfile', help='The full path to the seqfile.')
    parser.add_argument('--reference_file', help='The full path to the reference file for a cram.',
                        required=False, default=None)
    parser.add_argument('--threads', help='The number of threads to use.', type=int,
                        required=False, default=1)
    parser.add_argument('--skip_mapq', help='Skip plotting the mapq profile.', action='store_true')
    parser.add_argument('--skip_insert_sizes', help='Skip plotting the insert size  profile.',
                        action='store_true')
    params = parser.parse_args()

    if params.skip_mapq and params.skip_insert_sizes:
        raise RuntimeError('Cannot skip both mapq and insert size profile plotting.')

    seqfile = os.path.abspath(os.path.expanduser(params.seqfile))
    if not os.path.exists(params.seqfile):
        raise RuntimeError('Input seqfile (%s)does not exist.' % seqfile)

    reference_file=None
    if params.reference_file is not None:
        reference_file = os.path.abspath(os.path.expanduser(params.reference_file))
        if not os.path.exists(params.reference_file):
            raise RuntimeError('Input reference_file (%s)does not exist.' % reference_file)

    generate_mapplots(seqfile, params.threads, reference_file, params.skip_mapq,
                      params.skip_insert_sizes)

if __name__ == '__main__':
    main()
