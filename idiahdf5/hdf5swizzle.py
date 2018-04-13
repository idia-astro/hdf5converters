import logging
import numpy as np
import h5py

# TODO TODO TODO: move these to a different module called hdfswizzle.py
# Call that module from a script called hdf5swizzle
# Inside the script, autodetect XYZW -> ZYXW and use parallel for it -- call itself recursively with mpiexec; after handling all serial swizzles

class Dataset:
    pass # make this analogous to the converter dataset so that everything makes more sense

class Swizzler:
    PARALLEL = ('ZYXW',)
    
    def __init__(self, hdf5name):
        self.hdf5name = hdf5name
        
    def __enter__(self):
        self.hdf5 = h5py.File(self.hdf5name, "r+")
        return self
        
    def __exit__(self, e_type, e_value, e_traceback):
        self.hdf5.close()
    
    # TODO move to dataset
    def swizzle_axis_numeric(self, reversed_data_axes, axis_name):
        """Convert named swizzle axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, ZYXW -> (0, 3, 2, 1)
        """
        return tuple(reversed_data_axes.index(l) for l in reversed(axis_name))
    
    # TODO move to dataset
    def write_swizzled_dataset(self, axis_name):
        pass

    def swizzle(self, args):
        primary_group = self.hdf5["0"]
        data = primary_group["DATA"]
        
        swizzled_group = primary_group.require_group("SwizzledData")
        
        # TODO: this is good enough for now, but we should look for a more robust way to detect the order
        axes = "XYZW"[:len(data.shape)]
        reversed_data_axes = list(reversed(axes))
        
        for axis_name in args.swizzles:
            if not args.parallel or axis_name not in self.PARALLEL:
                

def swizzle(args):
    # do serial swizzles
    # do parallel swizzles
