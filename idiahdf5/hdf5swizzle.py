import logging
import numpy as np
import h5py
import subprocess

class Swizzler:
    def __init__(self, hdf5name):
        self.hdf5name = hdf5name
        self.hdfkwargs = {}
        
    def __enter__(self):
        self.hdf5 = h5py.File(self.hdf5name, "r+", **self.hdfkwargs)
        
        self.data = self.hdf5["0"]["DATA"]
        # TODO: this is good enough for now, but we should look for a more robust way to detect the order
        self.axes = "XYZW"[:len(self.data.shape)]
        self._reverse_axes = list(reversed(self.axes))
        
        return self
        
    def __exit__(self, e_type, e_value, e_traceback):
        self.hdf5.close()
    
    def swizzle_axis_numeric(self, reversed_data_axes, axis_name):
        """Convert named swizzle axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, ZYXW -> (0, 3, 2, 1)
        """
        return tuple(self._reverse_axes.index(l) for l in reversed(axis_name))
    
    def write_swizzled_dataset_serial(self, axis_name):
        if 'W' in axis_name and 'W' not in self.axes:
            axis_name = axis_name.replace('W', '')
            logging.warning("Trying to coerce swizzle axes to data axes. New swizzle axis: %s." % axis_name)
        
        if len(axis_name) != len(self.axes) or not all(d in self.axes for d in axis_name):
            logging.warning("Could not swizzle %s dataset to %s." % (self.axes, axis_name))
            return
        
        logging.info("Writing swizzled dataset %s..." % axis_name)
        
        axis = self.swizzle_axis_numeric(axis_name)
        swizzled = self.hdf5["0"].require_group("SwizzledData")
        
        swizzled.create_dataset(axis_name, data=np.transpose(self.data, axis), dtype=self.little_endian_dtype)
    
    def swizzle(self, args):
        self.parallel_todo = []
        
        for axis_name in args.swizzles:
            if args.parallel:
                parallel_params = SwizzlerParallel.implemented(self.data, axis_name)
                if parallel_params:
                    logging.info("Deferring %s swizzling to parallel implementation..." % axis_name)
                    self.parallel_todo.append(parallel_params) # leave until afterwards
                    continue # and skip serial swizzle
            
            logging.info("Swizzling %s..." % axis_name)
            self.write_swizzled_dataset_serial(axis_name) # otherwise do the serial swizzle


class SwizzlerParallel(Swizzler):
    @classmethod
    def implemented(cls, data, swizzled_axes):
        # only one case for now; may extend later
        if data.ndim == 4 and data.shape[0] > 1 and swizzled_axes == "ZYXW":
            return ('swizzle_ZYXW', data.shape[0]) # function to call and number of processes needed
        return None

    def __init__(self, hdf5name):
        super(SwizzlerParallel, self).__init__(hdf5name)

    def __enter__(self):
        from mpi4py import MPI
        self.hdf5kwargs = {'driver': 'mpio', 'comm': MPI.COMM_WORLD}
        super(SwizzlerParallel, self).__init__()
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
        

def swizzle(args):
    if args.quiet:
        logging.basicConfig(level=logging.CRITICAL)
    else:
        logging.basicConfig(level=logging.DEBUG)
    
    if args.funcname:
        # If we have a function name, this is being called recursively with MPI -- run the function and exit
        # We call each parallel swizzle individually, in case we want to execute them with different numbers of processes
        with SwizzlerParallel(args.filename) as swizzler:
            logging.debug("Process %d swizzling with function %s..." % (swizzler.rank, args.funcname))
            getattr(swizzler, args.funcname)()
            
    else:
        # Otherwise this is a normal call from the converter script or from the swizzle script
        
        # serial swizzles
        with Swizzler(args.filename) as swizzler:
            swizzler.swizzle(args)
            parallel_todo = swizzler.parallel_todo

        # parallel swizzles
        if args.parallel:
            for funcname, num_procs in parallel_todo:
                logging.info("Swizzling in parallel with %d processes and function %s..." % (num_procs, funcname))
                subprocess.run(['mpiexec', '-n', num_procs, 'hdf5swizzle', '--funcname', funcname, '--quiet', args.quiet])
