import argparse
from . import fits2hdf5

def hdf5_convert(filename, **kwargs):
    arguments = {
        "filename": filename,
        "statistics": ["XYZ", "XY", "Z"],
        "histograms": ["XYZ", "XY"],
        "percentiles": ["XYZ", "XY"],
        "swizzles": ["XYZW"],
        
        "quiet": True,
        "output_dir": None,
        "chunks": None,
    }
    
    arguments.update(kwargs)
    args = argparse.Namespace(**arguments)
    fits2hdf5.convert(args)
