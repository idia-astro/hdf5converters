import logging
import numpy as np
import h5py
import subprocess

# TODO simplify; merge this back into swizzler class
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
    def __init__(self, hdf5name):
        self.hdf5name = hdf5name
        
    def __enter__(self):
        self.hdf5 = h5py.File(self.hdf5name, "r+")
        self.data = self.hdf5["0"]["DATA"]
        return self
        
    def __exit__(self, e_type, e_value, e_traceback):
        self.hdf5.close()
    

    def swizzle(self, args):
        
        
        # TODO: this is good enough for now, but we should look for a more robust way to detect the order
        axes = "XYZW"[:len(data.shape)]
        reversed_data_axes = list(reversed(axes))
        
        for axis_name in args.swizzles:
            swizzled_group = self.hdf5["0"].require_group("SwizzledData")
            if not SwizzlerParallel.implemented(data, axis_name):
                pass # ??
                
                
                
class SwizzlerParallel(Swizzler):
    @classmethod
    def implemented(cls, data, swizzled_axes):
        # TODO: only one case for now; may extend later
        if data.ndim == 4 and data.shape[0] > 1 and swizzled_axes = "ZYXW":
            return True
        return False
    
    def __enter__(self):
        from mpi4py import MPI
        self.hdf5 = h5py.File(filename, 'r+', driver='mpio', comm=MPI.COMM_WORLD)
        self.data = self.hdf5["0"]["DATA"]
        self.rank = MPI.COMM_WORLD.rank
        return self
    
    def swizzle_ZYXW(self):
        W, Z, Y, X = self.data.shape
        
        swizzled_group = self.hdf5["0"].require_group("SwizzledData")
        swizzled_dataset = swizzled_group.create_dataset("ZYXW", (W, X, Y, Z), dtype=self.data.dtype)
        
        data_w = np.empty((Z, Y, X), dtype=self.data.dtype)
        
        w = self.rank
        
        # read one z at a time
        for z in range(Z):
            data_w[z] = self.data[w,z,:,:]
            
        # swizzle whole w by transposing
        swizzled_data_w = np.transpose(data_w, (2, 1, 0))
        
        # write one x at a time
        for x in range(X):
            swizzled_dataset[w, x, :, :] = swizzled_data_w[x, :, :]
            
    def swizzle(self):
        pass
                

def swizzle(args):
    # do serial swizzles
    # do parallel swizzles
