import logging
import numpy as np
import h5py

class Dataset:
    def __init__(self, hdu_group, data, axes):
        self.hdu_group = hdu_group
        self.axes = axes
        self._reverse_axes = list(reversed(axes))
        self.data = data
        self.little_endian_dtype = data.dtype.newbyteorder('L')

    def swizzle_axis_numeric(self, reversed_data_axes, axis_name):
        """Convert named swizzle axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, ZYXW -> (0, 3, 2, 1)
        """
        return tuple(self._reverse_axes.index(l) for l in reversed(axis_name))
    
    def write_swizzled_dataset(self, axis_name):
        if 'W' in axis_name and 'W' not in self.axes:
            axis_name = axis_name.replace('W', '')
            logging.warning("Trying to coerce swizzle axes to data axes. New swizzle axis: %s." % axis_name)
        
        if len(axis_name) != len(self.axes) or not all(d in self.axes for d in axis_name):
            logging.warning("Could not swizzle %s dataset to %s." % (self.axes, axis_name))
            return
        
        logging.info("Writing swizzled dataset %s..." % axis_name)
        
        axis = self.swizzle_axis_numeric(axis_name)
        swizzled = self.hdu_group.require_group("SwizzledData")
        
        swizzled.create_dataset(axis_name, data=np.transpose(self.data, axis), dtype=self.little_endian_dtype)


class Swizzler:
    PARALLEL = ('ZYXW',)
    
    def __init__(self, hdf5name):
        self.hdf5name = hdf5name
        
    def __enter__(self):
        self.hdf5 = h5py.File(self.hdf5name, "r+")
        return self
        
    def __exit__(self, e_type, e_value, e_traceback):
        self.hdf5.close()
    

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
