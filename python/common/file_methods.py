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
