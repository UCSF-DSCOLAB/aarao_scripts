import argparse
import gzip
import os
import re

from fastq_pair_iterator import *

fastq_re = re.compile(r'^(?P<prefix>(.*))_'
                      r'(?P<snum>(S[0-9]+))_'
                      r'(?P<lane>(L00[1-8]))_'
                      r'(?P<rtype>([IR][12]))_'
                      r'(?P<suffix>(001.fastq.gz))$')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq_dir', type=str, required=True, help='The path to a fastq folder.')
    parser.add_argument('--outdir', type=str, required=True, help='The path to an output directory'
                        '.')
    parser.add_argument('--remove_reads', type=str, required=True, help='A file containing '
                        'reads belonging to invalid samples (blacklist) with one read per line ('
                        'prefixed with @) and no header.')
    parser.add_argument('--library_name', type=str, required=False, help='The name of the library '
                        '(for output name) ', default=None)
    parser.add_argument('--update_verbosity', type=int, default=1_000_000)
    params = parser.parse_args()

    print('ACTION: Validating inputs')
    # params = parser.parse_args('--fastq_dir /krummellab/data1/staging/COMET_PBMC/fastqs/JYCVPBMC1/BCR1 --outdir . --remove_reads /krummellab/data1/staging/comet_PBMC_redacted/cr_outs/BCR_fastqs/MVIR1-POOL-SCB10/reads_to_remove.list.gz --library_name MVIR1-POOL-SCB10'.split())
    fastq_dir = os.path.abspath(params.fastq_dir)
    assert os.path.exists(fastq_dir)
    remove_reads_file = os.path.abspath(params.remove_reads)
    assert os.path.exists(remove_reads_file)

    outdir = os.path.abspath(params.outdir)

    if os.path.exists(outdir):
        outdir = os.path.join(outdir, '{}_fastqs'.format(params.library_name))
    else:
        assert os.path.exists(os.path.dirname(outdir))
    print('ACTION: Creating outdir {}'.format(outdir))
    os.mkdir(outdir)
    os.chdir(outdir)

    print('ACTION: Reading contents of {}'.format(remove_reads_file))
    assert is_gzipfile(remove_reads_file)
    remove_reads = {l.strip(b'\n') for l in gzip.open(remove_reads_file, 'rb')}

    fastq_files = [x for x in os.listdir(fastq_dir) if x.endswith('fastq.gz')]

    triplets = {}
    for ff in fastq_files:
        _ff = fastq_re.search(ff)
        assert _ff
        lane = _ff['lane'][-1].encode('utf-8')
        if lane not in triplets:
            triplets[lane] = {}
        triplets[lane][_ff['rtype']] = _ff
        assert is_gzipfile(os.path.join(fastq_dir, ff))
    assert all(len(x)==3 for x in triplets.values())

    remove_reads_per_lane = {}
    for r in remove_reads:
        l = r.split(b':')[-4]
        if l not in remove_reads_per_lane:
            remove_reads_per_lane[l] = set()
        remove_reads_per_lane[l].add(r)

    assert triplets.keys() == remove_reads_per_lane.keys()

    print('ACTION: Processing input fastq files....')
    total_reads = retain_reads = 0
    for lane in triplets:
        print('ACTION: Processing lane {} ....'.format(lane.decode('utf-8')))
        fp = paired_fastq(r1_file='{fastq_dir}/{prefix}_{snum}_{lane}_{rtype}_{suffix}'.format(fastq_dir=fastq_dir, **triplets[lane]['R1'].groupdict()),
                          r2_file='{fastq_dir}/{prefix}_{snum}_{lane}_{rtype}_{suffix}'.format(fastq_dir=fastq_dir, **triplets[lane]['R2'].groupdict()),
                          i1_file='{fastq_dir}/{prefix}_{snum}_{lane}_{rtype}_{suffix}'.format(fastq_dir=fastq_dir, **triplets[lane]['I1'].groupdict()))
        prefix2 = params.library_name or triplets[lane]['R1']['prefix']
        with gzip.open('{prefix2}_{snum}_{lane}_{rtype}_{suffix}'.format(prefix2=prefix2, **triplets[lane]['R1'].groupdict()), 'wb') as r1:
            with gzip.open('{prefix2}_{snum}_{lane}_{rtype}_{suffix}'.format(prefix2=prefix2, **triplets[lane]['R2'].groupdict()), 'wb') as r2:
                with gzip.open('{prefix2}_{snum}_{lane}_{rtype}_{suffix}'.format(prefix2=prefix2, **triplets[lane]['I1'].groupdict()), 'wb') as i1:
                    _total_reads = _retain_reads = 0
                    for reads in fp:
                        _total_reads += 1
                        rname = reads['R1'][0].split()[0]
                        l = r.split(b':')[-4]
                        if rname in remove_reads_per_lane[l]:
                            remove_reads_per_lane[l].remove(rname)
                            continue
                        r1.writelines(reads['R1'])
                        r2.writelines(reads['R2'])
                        i1.writelines(reads['I1'])
                        _retain_reads += 1
                        if _total_reads % params.update_verbosity == 0:
                            print('UPDATE: Processed {} read triplets.'.format(_total_reads))
        print('UPDATE: Processed {} read triplets.'.format(_total_reads))
        print('UPDATE: Retained {} read triplets.'.format(_retain_reads))
        total_reads += _total_reads
        retain_reads += _retain_reads
        fp.close()
    print('Process completed!')
    print('Processed a total of {} read triplets.'.format(total_reads))
    print('Retained a total of {} read triplets.'.format(retain_reads))

if __name__ == '__main__':
    main()