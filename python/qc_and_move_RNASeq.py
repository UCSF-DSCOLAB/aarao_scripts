import argparse
import glob
import os
import re
import sys

from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_directory')
    parser.add_argument('--action', choices=['qc', 'dryrun', 'final'], default='qc',
                        help="Just qc the data and print some stats, perform a dryrun, or actually move files.", 
                        required=False)
    parser.add_argument('--fastq_regex', default="^(?P<sample>(.*))_(?P<prefix>(S[0-9]{1,3}_L00[1-8]))_(?P<read>([IR][12]))_[0-9]{3}.fastq.gz", required=False)
    params = parser.parse_args()

    fastq_directory = os.path.abspath(params.fastq_directory)
    os.chdir(fastq_directory)

    qc = True
    dryrun = True
    if params.action == 'dryrun':
        qc = False
    if params.action == 'final':
        qc = False
        dryrun = False


    fastq_regex = re.compile(params.fastq_regex)
    if set(fastq_regex.groupindex.keys()) < {'sample', 'prefix', 'read'}:
        raise RuntimeError('fastq_regex must have 3 capture groups, `sample`, `prefix`, and `read`')
    files = glob.glob("*.fastq.gz")
    file_dict = defaultdict(dict)
    file_errors = 0
    for filename in files:
        if qc:
            print('processing: {}'.format(filename), end="\t")
        file_re = fastq_regex.search(filename)
        if not file_re:
            if qc:
                print('BAD_REGEX')
            file_errors += 1
            continue
        if file_re['prefix'] in file_dict[file_re['sample']]:
            if file_re['read'] in file_dict[file_re['sample']][file_re['prefix']]:
                print('DUPLICATE_READ_PREFIX')
                file_errors += 1
                continue
            else:
                if qc:
                    print('OK')

                file_dict[file_re['sample']][file_re['prefix']][file_re['read']]= filename
        else:
            if qc:
                print('OK')

            file_dict[file_re['sample']] = {
                file_re['prefix']: {
                    file_re['read']: filename
                }
            }

    if qc:
        print('processed all fastqs. Processing samples now...')
    sample_errors = 0
    for sample in file_dict:
        if qc:
            print('processing sample: {}'.format(sample))
        for prefix in file_dict[sample]:
            if qc:
                print('  processing prefix: {}{}'.format(sample, prefix), end="\t")
            if not {'R1', 'R2'}.issubset(file_dict[sample][prefix].keys()):
                sample_errors += 1
                if qc:
                    print('R1_OR_R2_MISSING')
            else:
                if qc:
                    print('OK')

    if qc or (file_errors + sample_errors):
        print('QC complete.\n'
              '    Found {} file errors\n'
              '    Found {} sample errors\n'.format(file_errors, sample_errors))
        sys.exit(file_errors > 0)

    for sample in file_dict:
        print('Creating folder: {}'.format(sample))
        if not dryrun:
            os.mkdir(sample)
        for prefix in file_dict[sample]:
            for filename in file_dict[sample][prefix]:
                print('Moving {} --> {}/{}'.format(filename, sample, filename))
                if not dryrun:
                    os.rename(filename, '{}/{}'.format(sample, filename))

if __name__ == '__main__':
    main()