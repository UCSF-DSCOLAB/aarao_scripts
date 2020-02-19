import os
import screed

class PairedScreedFastqs(object):
    def __init__(self, seqfile1, seqfile2):
        assert os.path.exists(seqfile1)
        assert os.path.exists(seqfile2)
        self.seqfile1 = screed.open(seqfile1)
        self.seqfile2 = screed.open(seqfile2)
    
    def __iter__(self): 
        return self

    def __next__(self):
        try:
            while True:
                return next(self.seqfile1.iter_fn), next(self.seqfile2.iter_fn)
        except StopIteration:
            print('Reached EOF.')
            return None, None

    def __del__(self):
        self.seqfile1.close()
        self.seqfile2.close()