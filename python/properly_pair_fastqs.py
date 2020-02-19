import argparse
import glob
import gzip
import os

from Bio import SeqIO


def is_gzipfile(filename):
    """
    Attempt to ascertain the gzip status of a file based on the "magic signatures" of the file.
    This was taken from the stack overflow post
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress
    :param str filename: A path to a file
    :return: True if the file appears to be gzipped else false
    :rtype: bool
    """
    assert os.path.exists(filename), 'Input {} does not point to a file.'.format(filename)
    with open(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == b'\x1f\x8b\x08':
            return True
        else:
            return False


def write_fastq(rname, r, f=None):
    if f is None:
        return
    print('@{}'.format(rname), file=f)
    print(r.seq, file=f)
    print('+', file=f)
    print(''.join(chr(x+33) for x in r.letter_annotations['phred_quality']), file=f)


def write_fastq_pair(rname, r1, r2, f1, f2):
    write_fastq(rname, r1, f1)
    write_fastq(rname, r2, f2)


def split_fastqs(if1, if2, of1, of2, ofu):
    u1 = {}
    u2 = {}
    
    r1 = r2 = None
    
    while True:
        r1 = next(if1, False)
        r2 = next(if2, False)
        if r1 and r2:
            if r1.id == r2.id:
                # Best case scenario
                write_fastq_pair(r1.id, r1, r2, of1, of2)
                continue
            if r1.id in u2:
                write_fastq_pair(r1.id, r1, u2.pop(r1.id), of1, of2)
                r1 = None
            if r2.id in u1:
                write_fastq_pair(r2.id, u1.pop(r2.id), r2, of1, of2)
                r2 = None
            if r1:
                u1[r1.id] = r1
            if r2:
                u2[r2.id] = r2
        elif not (r1 or r2):
            break
        elif not r2:
            u1[r1.id] = r1
            u1.update({r.id: r for r in if1})
            break
        else:  #elif not r1:
            u2[r2.id] = r2
            u2.update({r.id: r for r in if2})
            break

    for r in u1:
        if r in u2:
            write_fastq_pair(r, u1[r], u2.pop(r), of1, of2)
        else:
            write_fastq(r, u1[r], ofu)

    for r in u2:
        write_fastq(r, u2[r], ofu)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('F1', type=str)
    parser.add_argument('F2', type=str)
    parser.add_argument('--prefix', type=str, default=None, required=False)
    parser.add_argument('--R1', type=str, default=None, required=False)
    parser.add_argument('--R2', type=str, default=None, required=False)
    parser.add_argument('--RU', type=str, default=None, required=False)
    params = parser.parse_args()
    #params = parser.parse_args(['--samples', 'gyn_pr_myeloid.tsv'])

    F1 = os.path.abspath(params.F1)
    assert os.path.exists(F1)

    F2 = os.path.abspath(params.F2)
    assert os.path.exists(F2)

    R1, R2, RU = None, None, None
    if params.prefix:
        params.R1, params.R2, params.RU = [''.join(x) for x in zip([params.prefix]*3, ['_1.fq.gz', '_2.fq.gz', '_UP.fq.gz'])]
    else:
        assert any([params.R1, params.R2, params.RU]), "Need at least one of `prefix`, `R1`, `R2`, `RU`"
    
    F1 = (gzip.open if is_gzipfile(F1) else open)(F1, 'rt')
    F2 = (gzip.open if is_gzipfile(F2) else open)(F2, 'rt')
    
    R1 = gzip.open(os.path.abspath(params.R1), 'wt') if params.R1 else None
    R2 = gzip.open(os.path.abspath(params.R2), 'wt') if params.R2 else None
    RU = gzip.open(os.path.abspath(params.RU), 'wt') if params.RU else None

    try:
        split_fastqs(SeqIO.parse(F1, "fastq"), SeqIO.parse(F2, "fastq"), R1, R2, RU)
    finally:
        F1.close()
        F2.close()
        if R1 is not None:
            R1.close()
        if R2 is not None:
            R2.close()
        if RU is not None:
            RU.close()


if __name__ == '__main__':
    main()