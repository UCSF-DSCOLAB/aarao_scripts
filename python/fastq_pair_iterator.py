import gzip
import os

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


class paired_fastq(object):
    eof = False
    gzipped = False
    l_end = '\n'
    l_empty = ''
    r_start = '@'
    file_handles = None
    rstrip = False

    def __init__(self, r1_file=None, r2_file=None, i1_file=None, i2_file=None, gzip_mode='rb', rstrip=False):
        self.file_handles = {}
        self.gzipped = is_gzipfile(r1_file)
        if self.gzipped and 'b' in gzip_mode:
            self.l_end = b'\n'
            self.l_empty = b''
            self.r_start = b'@'
        self.rstrip = rstrip

        if is_gzipfile(r2_file) != self.gzipped:
            raise RuntimeError('Cannot mix and match compression status (all gzipped, or all not)')
        if i1_file and is_gzipfile(i1_file) != self.gzipped:
            raise RuntimeError('Cannot mix and match compression status (all gzipped, or all not)')
        if i2_file and is_gzipfile(i1_file) != self.gzipped:
            raise RuntimeError('Cannot mix and match compression status (all gzipped, or all not)')

        self.file_handles['R1'] = gzip.open(r1_file, gzip_mode) if self.gzipped else open(r1_file, 'r')
        self.file_handles['R2'] = gzip.open(r2_file, gzip_mode) if self.gzipped else open(r2_file, 'r')
        if i1_file:
            self.file_handles['I1'] = gzip.open(i1_file, gzip_mode) if self.gzipped else open(i1_file, 'r')
        if i2_file:
            self.file_handles['I2'] = gzip.open(i2_file, gzip_mode) if self.gzipped else open(i2_file, 'r')
    
    def close(self):
        self.file_handles['R1'].close()
        self.file_handles['R2'].close()
        if 'I1' in self.file_handles:
            self.file_handles['I1'].close()
        if 'I2' in self.file_handles:
            self.file_handles['I2'].close()

    def seek(self, offset, whence=0):
        if self.gzipped:
            raise RuntimeError("Cannot seek if the files are gzipped.")
        self.file_handles['R1'].seek(offset, whence)
        self.file_handles['R2'].seek(offset, whence)
        if 'I1' in self.file_handles:
            self.file_handles['I1'].seek(offset, whence)
        if 'I2' in self.file_handles:
            self.file_handles['I2'].seek(offset, whence)

    def next(self):
        if self.eof:
            raise StopIteration('EOF reached. Nothing to do.')
        record = {}
        temp = self.file_handles['R1'].readline()
        if temp == self.l_empty:
            self.eof = True
            raise StopIteration()

        assert temp.startswith(self.r_start)
        record['R1'] = [temp] + [self.file_handles['R1'].readline() for _ in range(3)]
        for r in self.file_handles:
            if r == 'R1':
                continue
            record[r] = [self.file_handles[r].readline() for _ in range(4)]
        if self.rstrip:
            record = {x: [_y.rstrip(self.l_end) for _y in y] for x, y in record.items()}
        return record
    
    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def __del__(self):
        self.close()