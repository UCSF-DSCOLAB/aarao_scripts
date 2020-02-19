import argparse
import json
import logging
import os
import re
import shutil
import sys
import yaml

from datetime import datetime as dt


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s: %(message)s')
logging.basicConfig(filemode='a', level=logging.DEBUG, format='%(levelname)s: %(message)s')

CTRL_REGEX = re.compile(r'CONTROL_(jurkat|UHR).Plate([0-9]{1,2}|Nova)')

compartments = {'live',
                'tcell',
                'treg',
                'myeloid',
                'tumor',
                'stroma',
                'epcam',
                'microglia',
                'total',
                'cd11bposhladrneg',
                'cd45pos',
                'cd45',
                'dpost',
                'tpost',
                'tpre',
                'cd45neg',
                'teff',
                'tcellp',
                'livep'}

IPI_REGEX = re.compile(r'^(?P<ipi>(IPI[A-Z]{2,4}[0-9]{3}))\.'
                       r'(?P<tn>(([TNRML]|NASH)[0-9]))\.'
                       r'rna\.'
                       r'(?P<stain>(' + '|'.join(compartments) + ')[2-9]*)$')


# Contains an approx tree of contents
expected_files = {
    'log': {
        '.trimmed.non_rrna.star.Log.final.out',
        '.trimmed.non_rrna.star.Log.out',
        'ht-pipes.log',
        },
    'metrics': {
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated.alignment_metrics',
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated_bamqc.pdf',
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated.flagstat',
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated.rnaseq_metrics',
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.duplication_metrics',
        '.trimmed.non_rrna.star.Aligned.toTranscriptome.out_bamqc.pdf',
        '.trimmed.rrna.sorted_bamqc.pdf',
        '.trimmed.rrna.sorted.flagstat',
        'fastp.html',
        'fastp.json',
        },
    'output': {
        '.rsem.genes.results',
        '.rsem.isoforms.results',
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated.cram',
        '.trimmed.non_rrna.star.Aligned.sortedByCoord.out.deduplicated.cram.crai',
        '.trimmed.non_rrna.star.Aligned.toTranscriptome.out.cram',
        '.trimmed.non_rrna.star.Chimeric.out.junction',
        '.trimmed.non_rrna.star.Unmapped.out.mate1.fastq.gz',
        '.trimmed.non_rrna.star.Unmapped.out.mate2.fastq.gz',
        '.trimmed.rrna.sorted.cram',
        '.trimmed.rrna.sorted.cram.crai',
        'fastq_md5sums.txt',
        },
    'results': {
        'gene_counts_table.tsv',
        'gene_exp_table.tsv',
        'gene_tpm_table.tsv',
        'isoform_tpm_table.tsv',
        'rnaseq_table.tsv',
        }
    }


def validate_plate_dir(plate_dir, old_sample, plate):
    """
    Validates the contents of `plate_dir`. Prints warnings for handle-able events and crashes on
    non-handleable ones.

    :param str plate_dir: The path to the plate directory
    :param dict(str, str|re) old_sample: A dict with the old sample 'name' and a 'regex' that
           matches to it
    :param str plate: the plate number as a string ('plate9', 'PlateNova', etc)
    """
    logger.info('Validating folder structure.')

    fail = 0

    main_folders = os.listdir(plate_dir)
    # Need to have a launch script and a yml
    if plate + '.rna_paired_align.sh' not in main_folders:
        fail += 1
        logging.error('{} does not contain a file named `{}`.'.format(
                      plate_dir, plate + '.rna_paired_align.sh'))
    if plate + '.rna_paired_align.yml' not in main_folders:
        fail += 1
        logging.error('{} does not contain a file named `{}`.'.format(
                      plate_dir, plate + '.rna_paired_align.yml'))

    # Get all samples from the yml
    with open(os.path.join(plate_dir, plate + '.rna_paired_align.yml')) as iff:
        yml = yaml.load(iff, Loader=yaml.FullLoader)
    samples = [yml[':samples'][i][':sample_name'] for i in range(len(yml[':samples']))]

    logger.debug('Detected %s samples:\nDEBUG\t%s', len(samples), '\nDEBUG:\t'.join(samples))
    # Check if the old sample regex matches multiple samples. This shouldn't happen and if it does,
    # I fucked up the regex
    samples_matching_old_regex = [s for s in samples if old_sample['regex'].search(s)]
    if len(samples_matching_old_regex) == 0:
        fail += 1
        logging.error('There are no samples in the yaml matching {}. Since we only modify the '
                      'yaml after successfully renaming all other files, this means something '
                      'is direly wrong, or your sample is not on this plate. Please continue '
                      'manually!'.format(old_sample['name']))
    elif len(samples_matching_old_regex) > 1:
        fail += 1
        logging.error('There are {} samples whose names start with {}. Safer to do this '
                      'manually'.format(len(samples_matching_old_regex), old_sample['name']))
    samples_matching_old_regex = samples_matching_old_regex[0]

    # Validate old and new names
    if old_sample['name'] not in samples:
        fail += 1
        logger.warning('%s does not contain a sample named `%s`.',
                       os.path.join(plate_dir, plate + '.rna_paired_align.yml'), old_sample['name'])

    # Now check folder structure.
    for _folder in ['metrics', 'output', 'log']:
        if _folder not in main_folders:
            fail += 1
            logging.error('%s does not contain a folder named `%s`.', plate_dir, _folder)
            continue

        folder = os.path.join(plate_dir, _folder)
        if old_sample['name'] not in os.listdir(folder):
            fail += 1
            logger.warning('%s does not contain a folder named `%s`.', folder, old_sample['name'])
            continue
        else: # old_sample['name'] exists and new_sample['name'] doesn't.
            folder = os.path.join(plate_dir, folder, old_sample['name'])

        files, old_files= set(), set()
        for f in os.listdir(folder):
            if f.startswith('.'):
                continue

            if old_sample['regex'].search(f):
                old_files.add(f)
                f = old_sample['regex'].sub('', f)
            else:
                old_files.add(f)
            files.add(f)

        logger.debug('Old files in %s:\nDEBUG\t%s', folder, '\nDEBUG:\t'.join(old_files))


        if files != expected_files[_folder]:
            fail += 1
            logger.error('The contents of %s did not match the expected contents', folder)
            logger.error('Expected contents: %s', ','.join(expected_files[_folder]))
            logger.error('Observed contents: %s', ','.join(files))
            
    if 'results' not in main_folders:
        logger.info('%s did not contain a folder named `results` but that\'s ok', plate_dir)
    else:
        files = {x for x in os.listdir(os.path.join(plate_dir, 'results')) if not x.startswith('.')}
        if files != expected_files['results']:
            fail += 1
            logger.error('The contents of %s did not match the expected contents',
                          os.path.join(plate_dir, results))
            logger.error('Expected contents: %s', ','.join(expected_files['results']))
            logger.error('Observed contents: %s', ','.join(files))
            
    if 'scripts' not in main_folders:
        fail += 1
        logging.error('{} does not contain a folder named `scripts`.'.format(plate_dir))
    else:
        files = {x for x in os.listdir(os.path.join(plate_dir, 'scripts'))
                 if old_sample['regex'].search(x)}
        if len(files) != 1:
            assert len(files) == 0, 'WTF1'
            fail += 1
            logger.warning('{} does not contain a single file matching `{}`.'.format(
                          os.path.join(plate_dir, 'scripts'), old_sample['regex']))

    logger.warning('Detected %s warnings!!!', fail)
    return fail


def main():
    """
    This program provides a utility to validate a single RNASeq sample on an IPI/ImmunoX RNASeq plate.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--plate_dir', type=str, required=True)
    parser.add_argument('--plate', type=str, required=True, help='E.g.: plateNova, plate11, etc.')
    parser.add_argument('--sample_name', type=str, required=True)
    params = parser.parse_args()

    start = dt.now()
    if not os.path.exists(params.plate_dir):
        logger.info(1)
        raise RuntimeError('Non-existant `plate_dir`: {}'.format(params.plate_dir))
    else:
        plate_dir = os.path.abspath(params.plate_dir)

    # In the regexes, convert periods to [.] so they are interpreted literally. Add a lookahead so
    # the replace can be used to differentiate between `IPIREF999.T1.rna.tcell` and
    # `IPIREF999.T1.rna.tcell2`
    old_sample = {
        'name': params.sample_name,
        'regex': re.compile(params.sample_name.replace('.', '[.]') + '(?![p0-9])')
    }

    logger.info(dt.strftime(start,
                'Processess started at %H:%M:%S on %B %d, %Y'))
    logger.info('Running with command:\n\t\t%s', ' '.join(sys.argv))
    logger.info('\n\t'.join([
        'Parameters:',
        'plate_dir: %s' % params.plate_dir,
        'plate: %s' % params.plate,
        'sample_name: %s' % params.sample_name
        ]))

    fail = validate_plate_dir(plate_dir, old_sample, params.plate)
    logger.info('Process took %s' % str(dt.now()-start).split('.', 2)[0])
    logger.info(fail)

if __name__ == '__main__':
    main()