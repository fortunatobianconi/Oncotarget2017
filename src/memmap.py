from os import path
import tempfile
import gc
from joblib import Parallel, delayed, load, dump

def memmap(a):
    '''use memory mapping to prevent memory allocation for each worker'''
    tmp_dir = tempfile.mkdtemp()
    mmap_fn = path.join(tmp_dir, 'a.mmap')
    print ('mmap file:', mmap_fn)
    _ = dump(a, mmap_fn)        # dump
    a_mmap = load(mmap_fn, 'r+') # load
    del a
    gc.collect()
    return a_mmap