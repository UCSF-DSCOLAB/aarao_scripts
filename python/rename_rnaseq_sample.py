import argparse
import json
import logging
import os
import re
import pysam
import shutil
import sys
import yaml

from datetime import datetime as dt


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s: %(message)s')

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


def validate_plate_dir(plate_dir, old_sample, new_sample, plate, resume):
    """
    Validates the contents of `plate_dir`. Prints warnings for handle-able events and crashes on
    non-handleable ones.

    :param str plate_dir: The path to the plate directory
    :param dict(str, str|re) old_sample: A dict with the old sample 'name' and a 'regex' that
           matches to it
    :param dict(str, str|re) new_sample: A dict with the new sample 'name' and a 'regex' that
           matches to it
    :param str plate: the plate number as a string ('plate9', 'PlateNova', etc)
    :param bool resume: Is this a resumed run?
    """
    logger.info('Validating folder structure.')

    fail = 0  # To tag events that don't require an automatic raise RuntimeError

    _new_sample_name = IPI_REGEX.search(new_sample['name']) or CTRL_REGEX.search(new_sample['name'])
    if not _new_sample_name:
        raise RuntimeError('The new sample name is invalid.')

    main_folders = os.listdir(plate_dir)
    # Need to have a launch script and a yml
    if plate + '.rna_paired_align.sh' not in main_folders:
        raise RuntimeError('{} does not contain a file named `{}`.'.format(
            plate_dir, plate + '.rna_paired_align.sh'))
    if plate + '.rna_paired_align.yml' not in main_folders:
        raise RuntimeError('{} does not contain a file named `{}`.'.format(
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
        raise RuntimeError('There are no samples in the yaml matching {}. Since we only modify the '
                           'yaml after successfully renaming all other files, this means something '
                           'is direly wrong, or your sample is not on this plate. Please continue '
                           'manually!'.format(old_sample['name']))
    elif len(samples_matching_old_regex) > 1:
        raise RuntimeError('There are {} samples whose names start with {}. Safer to do this '
                           'manually'.format(len(samples_matching_old_regex), old_sample['name']))
    samples_matching_old_regex = samples_matching_old_regex[0]
    # Check if the old sample regex matches the fastq files because if they do, we will definitely
    # fuck up fastq_md5sums.txt, ht-pipes.log, etc
    for fq_file in yml[':samples'][samples.index(samples_matching_old_regex)][':sample_files']:
        if old_sample['regex'].search(fq_file):
            raise RuntimeError('The old sample regex matches an input fastq file (%s) so you '
                               'really should do this renaming thing manually.' % fq_file)

    # Validate old and new names
    if old_sample['name'] not in samples:
        logger.warning('%s does not contain a sample named `%s`.',
                       os.path.join(plate_dir, plate + '.rna_paired_align.yml'), old_sample['name'])
        fail += 1
        if new_sample['name'] in samples:
            logger.warning('However %s contains a folder named `%s`',
                            os.path.join(plate_dir, plate + '.rna_paired_align.yml'),
                            new_sample['name'])
        else:
            raise RuntimeError('{} does not contain a folder named `{}` either.'.format(
                os.path.join(plate_dir, plate + '.rna_paired_align.yml'), new_sample['name']))
    elif new_sample['name'] in samples:  # Implies old_sample['name'] exists too.
        raise RuntimeError('{} already contains both, samples named `{}`, and `{}`.'.format(
            os.path.join(plate_dir, plate + '.rna_paired_align.yml'), old_sample['name'],
            new_sample['name']))

    # Now check folder structure.
    for _folder in ['metrics', 'output', 'log']:
        if _folder not in main_folders:
            raise RuntimeError('%s does not contain a folder named `%s`.', plate_dir, _folder)

        folder = os.path.join(plate_dir, _folder)
        if old_sample['name'] not in os.listdir(folder):
            fail += 1
            logger.warning('%s does not contain a folder named `%s`.', folder, old_sample['name'])
            if new_sample['name'] in os.listdir(folder):
                logger.warning('However %s contains a folder named `%s`.',
                                folder, new_sample['name'])
                folder = os.path.join(plate_dir, folder, new_sample['name'])
            else:
                raise RuntimeError('{} does not contain a folder named `{}` either.'.format(
                    folder, new_sample['name']))
        elif new_sample['name'] in os.listdir(folder):  # old_sample['name'] exists.
            raise RuntimeError('{} already contains both, samples named `{}`, and `{}`.'.format(
                folder, old_sample['name'], new_sample['name']))
        else: # old_sample['name'] exists and new_sample['name'] doesn't.
            folder = os.path.join(plate_dir, folder, old_sample['name'])

        files, old_files, new_files = set(), set(), set()
        for f in os.listdir(folder):
            if f.startswith('.'):
                continue

            if old_sample['regex'].search(f):
                old_files.add(f)
                f = old_sample['regex'].sub('', f)
            elif new_sample['regex'].search(f):
                new_files.add(f)
                f = new_sample['regex'].sub('', f)
            else:
                old_files.add(f)
            files.add(f)

        logger.debug('Old files in %s:\nDEBUG\t%s', folder, '\nDEBUG:\t'.join(old_files))
        logger.debug('Updated files in %s:\nDEBUG\t%s', folder, '\nDEBUG:\t'.join(new_files))


        if files != expected_files[_folder]:
            logger.error('The contents of %s did not match the expected contents', folder)
            logger.error('Expected contents: %s', ','.join(expected_files[_folder]))
            logger.error('Observed contents: %s', ','.join(files))
            assert False
        if len(new_files) > 0:
            fail += 1
            logger.warning('%s had %s files starting with %s', folder, len(new_files),
                            new_sample['name'])

    if 'results' not in main_folders:
        logger.info('%s did not contain a folder named `results` but that\'s ok', plate_dir)
    else:
        files = {x for x in os.listdir(os.path.join(plate_dir, 'results')) if not x.startswith('.')}
        if files != expected_files['results']:
            logger.error('The contents of %s did not match the expected contents',
                          os.path.join(plate_dir, 'results'))
            logger.error('Expected contents: %s', ','.join(expected_files['results']))
            logger.error('Observed contents: %s', ','.join(files))
            assert False

    if 'scripts' not in main_folders:
        raise RuntimeError('{} does not contain a folder named `scripts`.'.format(plate_dir))
    else:
        files = {x for x in os.listdir(os.path.join(plate_dir, 'scripts'))
                 if old_sample['regex'].search(x)}
        if len(files) != 1:
            assert len(files) == 0, 'WTF1'
            fail += 1
            logger.warning('{} does not contain a single file matching `{}`.'.format(
                          os.path.join(plate_dir, 'scripts'), old_sample['regex']))
            files = {x for x in os.listdir(os.path.join(plate_dir, 'scripts'))
                     if new_sample['regex'].search(x)}
            if len(files) != 1:
                assert len(files) == 0, 'WTF2'
                logger.error('{} does not contain a single file matching `{}` either.'.format(
                              os.path.join(plate_dir, 'scripts'), new_sample['regex']))
            else:
                logger.warning('However, {} contains a single file matching `{}`.'.format(
                              os.path.join(plate_dir, 'scripts'), new_sample['regex']))

    if fail > 0:
        logger.warning('Detected %s warnings!!!', fail)
        if resume:
            logger.warning('Folder structure had some issues but we\'re going to continue since '
                            'the `--resume flag was provided`!')
        else:
            raise RuntimeError('Encountered an error with the folder structure on a fresh run. '
                               'If you are attempting to resume a previously killed run, use the '
                               '`--resume` flag.')
    else:
        logger.info('Folder structure looks good!')


def modify_text_file(old_file, old_sample, new_sample, dry):
    """
    Modifies the contents of a text file to change the old sample name(using a regex) to the new one
    and writes it to a temp file. Then renames the temp file to the new file name before deleting
    the old file.

    :param str old_file: The full path to the old file
    :param dict(str, str|re) old_sample: A dict with the old sample 'name' and a 'regex' that
           matches to it
    :param dict(str, str|re) new_sample: A dict with the new sample 'name' and a 'regex' that
           matches to it
    :param bool dry: Is this a dry run?
    """
    old_file = os.path.abspath(old_file)
    contents = [line for line in open(old_file)]

    for i, _ in enumerate(contents):
        if old_sample['regex'].search(contents[i]):
            logger.info('Modifying line %s' % (i+1))
            logger.debug('OLD: %s', contents[i])
            contents[i] = old_sample['regex'].sub(new_sample['name'], contents[i])
            logger.debug('NEW: %s', contents[i])

    new_file = os.path.join(os.path.dirname(old_file),
                                            old_sample['regex'].sub(new_sample['name'],
                                                                    os.path.basename(old_file)))
    temp_file = os.path.join(os.path.dirname(new_file), '.' + os.path.basename(new_file))
    logger.debug('Writing %s', os.path.basename(new_file))
    if not dry:
        logger.debug('Writing temp %s', os.path.basename(temp_file))
        with open(temp_file, 'w') as off:
            print(''.join(contents), file=off, end='')
        logger.debug('Actually writing %s', os.path.basename(new_file))
        shutil.move(temp_file, new_file)
        if not old_file == new_file:
            # Names without the sample prefix don't get renamed.
            os.remove(old_file)


def modify_cram(old_file, old_sample, new_sample, ref_fasta, dry):
    """
    Modifies the header of a cram file to change the old sample name(using a regex) to the new one
    and writes it to a temp file. Then renames the temp file to the new file name before deleting
    the old file and old index (if present). The new cram is indexed if required.

    :param str old_file: The full path to the old file
    :param dict(str, str|re) old_sample: A dict with the old sample 'name' and a 'regex' that
           matches to it
    :param dict(str, str|re) new_sample: A dict with the new sample 'name' and a 'regex' that
           matches to it
    :param str ref_fasta: Path to the reference 'genome', 'transcriptome', or  'rrna' fasta for the
           cram file
    :param bool dry: Is this a dry run?
    """
    old_file = os.path.abspath(old_file)
    changes = 0
    with pysam.AlignmentFile(old_file, reference_filename=ref_fasta) as old_sample_cram:
        header = old_sample_cram.header.as_dict()
        index = old_sample_cram.has_index()
        if 'PG' in header:
            for i in range(len(header['PG'])):
                logger.info('Modifying item %s in the PG tag' % (i+1))
                logger.debug('OLD: %s', header['PG'][i]['CL'])
                header['PG'][i]['CL'] = old_sample['regex'].sub(new_sample['name'],
                                                                header['PG'][i]['CL'])
                logger.debug('NEW: %s', header['PG'][i]['CL'])
                changes += 1
        if 'CO' in header:
            logger.info('Modifying the CO tag')
            logger.debug('OLD: %s', header['CO'])
            for i in range(len(header['CO'])):
                header['CO'][i] = old_sample['regex'].sub(new_sample['name'], header['CO'][i])
                changes += 1
            logger.debug('NEW: %s', header['CO'])

        new_file = os.path.join(os.path.dirname(old_file),
                                old_sample['regex'].sub(new_sample['name'],
                                                        os.path.basename(old_file)))
        temp_file = os.path.join(os.path.dirname(new_file), '.' + os.path.basename(new_file))

        logger.info('Creating new cram file %s', new_file)
        if not dry:
            if changes > 0 :
                logger.info('Creating temp %s', temp_file)
                with pysam.AlignmentFile(temp_file, 'wc', reference_filename=ref_fasta,
                                         header=header) as new_sample_cram:
                    for s in old_sample_cram:
                        new_sample_cram.write(s)
                logger.debug('Actually writing %s', os.path.basename(new_file))
                shutil.move(temp_file, new_file)
                os.remove(old_file)
            else:
                logger.info('No changes detected. Creation is via renaming.')
                shutil.move(old_file, new_file)

    if index:
        logger.info('Creating new cram index file for %s', new_file)
        if not dry:
            pysam.index(new_file, new_file + '.crai')
            os.remove(old_file + '.crai')
    else:
        logger.info('Cram file %s does not need an index', new_file)


def process_sample(plate_dir, old_sample, new_sample, plate, fastas, dry):
    """
    Processes all files in the plate directory associated with the sample including modifying the
    file contents to change the old sample name to the new one, and renaming the file (if required).

    Rename operations are atomic and files that need to be copied are first written to a temporary
    file and then atomically renamed to the new name before the old file is removed. This way we
    preserve the original file till the new one is created and there is never a state where the new
    one is incomplete and the old one is removed/modified.

    :param str plate_dir: The path to the plate directory
    :param dict(str, str|re) old_sample: A dict with the old sample 'name' and a 'regex' that
           matches to it
    :param dict(str, str|re) new_sample: A dict with the new sample 'name' and a 'regex' that
           matches to it
    :param str plate: the plate number as a string ('plate9', 'PlateNova', etc)
    :param dict(str, str) fastas: A map to the 'genome', 'transcriptome', and 'rrna' fasta paths
    :param bool dry: Is this a dry run?
    """
    logger.info('Entering folder %s' % plate_dir)
    for folder in expected_files:
        _folder = os.path.join(plate_dir, folder)
        _folder_new = None
        if folder == 'results':
            if not os.path.exists(_folder):
                # This is ok
                continue
        else:
            if old_sample['name'] in os.listdir(_folder):
                _folder_new = os.path.join(_folder, new_sample['name'])
                _folder = os.path.join(_folder, old_sample['name'])
            elif new_sample['name'] in os.listdir(_folder):
                _folder = os.path.join(_folder, new_sample['name'])
            else:
                assert False, 'Structure validated but neither new nor old samples in this folder?'
        logger.info('Entering folder: %s' % _folder)

        for f in os.listdir(_folder):
            if f.startswith('.'):
                continue
            logger.info('Processing file: %s' % f)
            if new_sample['regex'].search(f):
                # We don't name a file till we fill its contents so this can safely be assumed
                # to have been handled
                logger.info('Detected file is already processed (based on file name). Skipping...')
                continue

            if f.endswith('cram'):
                _f = ('genome' if 'sortedByCoord' in f
                      else 'transcriptome' if 'toTranscriptome' in f
                      else 'rrna')
                modify_cram(os.path.join(_folder, f), old_sample, new_sample, fastas[_f], dry)
            else:
                if f.endswith('crai'):
                    # Nothing to do. Handled in the cram case.
                    continue
                if f.endswith(('fastq.gz', 'pdf')):
                    # Just rename
                    f_new = old_sample['regex'].sub(new_sample['name'], f)
                    logger.info('Renaming to %s', f_new)
                    if not dry:
                        shutil.move(os.path.join(_folder, f), os.path.join(_folder, f_new))
                else:
                    modify_text_file(os.path.join(_folder, f), old_sample, new_sample, dry)
        if _folder_new is not None:
            logger.info('Renaming the entire folder to %s', _folder_new)
            if not dry:
                shutil.move(_folder, _folder_new)
        else:
            logger.info('This folder does not require to be renamed.')
        logger.info('Leaving folder: %s', _folder)

    files = {x for x in os.listdir(os.path.join(plate_dir, 'scripts'))
             if old_sample['regex'].search(x)}
    for f in files:
        logger.info('Entering folder: scripts')
        logger.info('Processing file %s', f)
        modify_text_file(os.path.join(plate_dir, 'scripts', f), old_sample, new_sample, dry)
        logger.info('Leaving folder: scripts')

    logger.info('Finally... Processing file %s', plate + '.rna_paired_align.yml')
    modify_text_file(os.path.join(plate_dir, plate + '.rna_paired_align.yml'), old_sample,
                     new_sample, dry)


def main():
    """
    This program provides a utility to rename a single RNASeq sample on an IPI/ImmunoX RNASeq plate.
    We validate the existence of the file in the provided `plate_dir` and the rename files and
    modify file internals to make it seem like the file was named `new_sample_name` all along. The
    log of the renaming including every action taken is stored in a unique file within `plate_dir`.

    There is an option to do a "Dry Run" (`--dry`) so you can personally see what changes will be
    made during the process. All naming operations are atomic so this script is extremely robust to
    restarting.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--plate_dir', type=str, required=True)
    parser.add_argument('--plate', type=str, required=True, help='E.g.: plateNova, plate11, etc.')
    parser.add_argument('--old_sample_name', type=str, required=True)
    parser.add_argument('--new_sample_name', type=str, required=True)
    parser.add_argument('--resume', action='store_true', help='Resume a previously cancelled run?')
    parser.add_argument('--dry', action='store_true', help='Dry run. Don\'t rename anything.')
    parser.add_argument('--genome_fa', help='Path to the genome fasta.',
                        default='/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa')
    parser.add_argument('--rrna_fa', help='Path to the rrna fasta.',
                        default='/krummellab/data1/ipi/data/refs/bwa/human.rrna.fa')
    parser.add_argument('--transcriptome_fa', help='Path to the transcriptome fasta.',
                        default='/krummellab/data1/ipi/data/refs/rsem/hg38_pure_rsem/hg38.transcripts.fa')
    parser.add_argument('--log_to_console', action='store_true', help='Log to STDOUT?.  All '
                        'events are automatically logged to a file in the plate directory.')
    parser.add_argument('--log_level', type=str, choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help='The level of logging above which messages should be '
                        'printed to console.')
    params = parser.parse_args()

    start = dt.now()
    if not os.path.exists(params.plate_dir):
        raise RuntimeError('Non-existant `plate_dir`: {}'.format(params.plate_dir))
    else:
        plate_dir = os.path.abspath(params.plate_dir)

    fastas = {
        'genome': '',
        'rrna': '',
        'transcriptome': ''
    }
    for f in fastas:
        _f = getattr(params, '%s_fa' % f)
        if not os.path.exists(_f):
            raise RuntimeError('Non-existant `{}_fa`: {}'.format(f, _f))
        else:
            fastas[f] = os.path.abspath(_f)

    # In the regexes, convert periods to [.] so they are interpreted literally. Add a lookahead so
    # the replace can be used to differentiate between `IPIREF999.T1.rna.tcell` and
    # `IPIREF999.T1.rna.tcell2`
    old_sample = {
        'name': params.old_sample_name,
        'regex': re.compile(params.old_sample_name.replace('.', '[.]') + '(?![a-zA-Z0-9])')
    }
    new_sample = {
        'name': params.new_sample_name,
        'regex': re.compile(params.new_sample_name.replace('.', '[.]') + '(?![a-zA-Z0-9])')
    }

    # Ensure the regexes don't match to each other
    if old_sample['name'] == new_sample['name']:
        raise RuntimeError('Old and new samples cannot be the same')
    if old_sample['regex'].search(new_sample['name']):
        raise RuntimeError('The old sample regex matches the new sample name. Not advisable')
    if new_sample['regex'].search(old_sample['name']):
        raise RuntimeError('The new sample regex matches the old sample name. Not advisable')

    # Create logfiles only if the input arguments are valid
    logfile = os.path.join(plate_dir,
                           dt.strftime(dt.now(), '%y_%m_%d__%H_%M_%S_rename_samples.log'))
    if params.dry:
        logfile = re.sub('.log$', '_DRY.log', logfile)
    print('Logging to file: %s' % logfile)
    lfh = logging.FileHandler(logfile)
    lfh.setFormatter(formatter)
    lfh.setLevel(logging.DEBUG)
    logger.addHandler(lfh)

    if params.log_to_console:
        print('Additionally logging to STDOUT.')
        console = logging.StreamHandler()
        console.setLevel(getattr(logging, params.log_level))
        console.setFormatter(formatter)
        logger.addHandler(console)

    logger.info(dt.strftime(start,
                'Processess started at %H:%M:%S on %B %d, %Y'))
    logger.info('Running with command:\n\t\t%s', ' '.join(sys.argv))
    logger.info('\n\t'.join([
        'Parameters:',
        'plate_dir: %s' % params.plate_dir,
        'plate: %s' % params.plate,
        'old_sample_name: %s' % params.old_sample_name,
        'new_sample_name: %s' % params.new_sample_name,
        'resume: %s' % params.resume,
        'dry: %s' % params.dry,
        'genome_fa: %s' % params.genome_fa,
        'rrna_fa: %s' % params.rrna_fa,
        'transcriptome_fa: %s' % params.transcriptome_fa,
        'log_to_console: %s' % params.log_to_console,
        'log_level: %s' % params.log_level
        ]))

    lock_file = os.path.join(plate_dir, '.lock')
    lockfile_contents = [old_sample['name'], 
                         new_sample['name'], 
                         dt.strftime(start, '%H:%M:%S  %B %d, %Y')]

    if os.path.exists(lock_file):
        raise RuntimeError('Found the lock file at %s. This usually means you have a running rename '
                           'script on this directory so you shouldn\'t be starting a new one till '
                           'that one is done. If you know you don\'t have a running operation, you '
                           'can delete the lock file and continue (at your own risk). The lock '
                           'file contents will tell you what sample is being renamed.' % lock_file)
    else:
        logger.info('Creating a lock file at %s to prevent other processes from running '
                     'simultaneously' % lock_file)
        with open(lock_file, 'w') as iff:
            print('\n'.join(lockfile_contents), file=iff)
    try:
        validate_plate_dir(plate_dir, old_sample, new_sample, params.plate, params.resume)
        process_sample(plate_dir, old_sample, new_sample, params.plate, fastas, params.dry)
    except:
        logger.error(dt.strftime(dt.now(),
                     'Process failed at %H:%M:%S on %B %d, %Y'))
        raise
    else:
        logger.info(dt.strftime(dt.now(),
                    'Process completed successfully at %H:%M:%S on %B %d, %Y'))
    finally:
        logger.info('Process took %s' % str(dt.now()-start).split('.', 2)[0])
        contents = [l.strip() for l in open(lock_file)]
        if contents != lockfile_contents:
            logger.error('Lockfile contents do not match the expected ones. This can suggest '
                          'this run has been compromised. Not removing the existing lockfile '
                          'but ¯\\_(-_-)_/¯')
            logger.error('Contents: %s', ', '.join(contents))
            logger.error('Expected: %s', ', '.join(lockfile_contents))
            raise RuntimeError('Lock File Content Mismatch')
        os.remove(lock_file)


if __name__ == '__main__':
    main()
