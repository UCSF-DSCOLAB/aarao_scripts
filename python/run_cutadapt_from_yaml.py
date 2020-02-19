from __future__ import print_function

import argparse
import glob
import os
import subprocess
import yaml

def get_cutadapt_stats(input_yaml_file, output_folder, kit):
    assert kit in ['nextera', 'truseq']
    with open(input_yaml_file) as iff:
        input_yaml = yaml.load(iff)
    plate = os.path.basename(input_yaml_file).split('.')[0]
    _output_folder = os.path.join(output_folder, plate)
    if not os.path.exists(_output_folder):
        os.mkdir(_output_folder, 0777)
    print('Processing %s' % plate)
    _plate = plate[0].upper() + plate[1:]
    if os.path.exists(os.path.join('/cbc2/data2/samadb/IPI/fastqs', _plate + 'Rnaseq')):
        fastq_dir = os.path.join('/cbc2/data2/samadb/IPI/fastqs', _plate + 'Rnaseq')
    elif os.path.exists(os.path.join('/cbc2/data1/data/collaborations/Krummel/', _plate + 'Rnaseq')):
        fastq_dir = os.path.join('/cbc2/data1/data/collaborations/Krummel/', _plate + 'Rnaseq')
    else:
        print('Cannot identify fastq dir for %s' % plate)
        return None
    print('Identified fastq dir for %s as %s' % (plate, fastq_dir))

    samples = {x[':sample_name']: (os.path.basename(x[':replicates'][0][':inputs'][0][':fq1']), os.path.basename(x[':replicates'][0][':inputs'][0][':fq2'])) for x in input_yaml[':samples']}
    
    plate_files = {os.path.basename(x): x for x in glob.glob(os.path.join(fastq_dir,'*','*.gz'))}

    for sample in samples:
        if not (samples[sample][0] in plate_files and samples[sample][1] in plate_files):
            print('ERROR:', sample, samples[sample][0] in plate_files, samples[sample][1] in plate_files, sep='\t')
            continue
        call = [
            '/home/arrao/.local/bin/cutadapt',
            '-j', os.environ['PBS_NUM_PPN'] if 'PBS_NUM_PPN' in os.environ else '1',
            '-o', '/dev/null',
            '-p', '/dev/null',
            plate_files[samples[sample][0]],
            plate_files[samples[sample][1]]
        ]
        if kit == 'truseq':
            call.extend([
                '-a', 'AGATCGGAAGAG',
                '-A', 'AGATCGGAAGAG'
            ])
        else:
            call.extend([
                '-a', 'CTGTCTCTTATACACATCT',
                '-A', 'CTGTCTCTTATACACATCT'
            ])

        base = os.path.splitext(os.path.basename(input_yaml_file))[0]
        with open(os.path.join(_output_folder, base + '_cutadapt_stats.txt'), 'w') as off:
            try:
                subprocess.check_call(call, stdout=off)
            except subprocess.CalledProcessError as e:
                print('ERROR: Failed on %s:%s' % (plate, sample))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_yaml', help='Input yaml.', type=str)
    parser.add_argument('--output_folder', help='Output folder.', type=str,
                        default='/data/shared/krummellab/ipi/cutadapt_stats', required=False)
    parser.add_argument('--kit', help='Illumina Kit used.', choices=['nextera', 'truseq'],
                        type=str, default='truseq')
    params, extras = parser.parse_known_args()


    assert os.path.exists(params.input_yaml), 'Input yaml does not exist.'
    assert os.path.exists(params.output_folder), 'Output folder does not exist.'

    get_cutadapt_stats(params.input_yaml, params.output_folder, params.kit)

if __name__ == '__main__':
    main()