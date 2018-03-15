import argparse
from .fits2hdf5 import convert as f2h_convert

def convert(filename, **kwargs):
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
    f2h_convert(args)
