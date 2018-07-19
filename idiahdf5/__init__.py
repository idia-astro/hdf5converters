import argparse
from .fits2hdf5 import convert, DEFAULTS

def hdf5_convert(filename, **kwargs):
    arguments = {
        "filename": filename,
        "statistics": DEFAULTS.statistics,
        "histograms": DEFAULTS.histograms,
        "percentiles": DEFAULTS.percentiles,
        "swizzles": DEFAULTS.swizzles,
        
        "quiet": True,
        "output_dir": None,
        "chunks": None,
    }
    
    arguments.update(kwargs)
    args = argparse.Namespace(**arguments)
    convert(args)
