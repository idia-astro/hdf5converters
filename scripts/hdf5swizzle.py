#!/usr/bin/env python3

import h5py
import numpy as np
from mpi4py import MPI
import sys

comm = MPI.COMM_WORLD
rank = comm.rank
nprocs = comm.size

def swizzle(filename, funcname): # for this prototype assume XYZW -> ZYXW
    with h5py.File(filename, 'r+', driver='mpio', comm=MPI.COMM_WORLD) as f:
        data = f["0"]["DATA"]
        W, Z, Y, X = data.shape
        
        swizzled_group = f["0"].require_group("SwizzledData")
        swizzled_dataset = swizzled_group.create_dataset("ZYXW", (W, X, Y, Z), dtype=data.dtype)
        
        data_w = np.empty((Z, Y, X), dtype=data.dtype)
        
        # read one z at a time
        for z in range(Z):
            data_w[z] = data[w,z,:,:]
            
        # swizzle whole w by transposing
        swizzled_data_w = np.transpose(data_w, (2, 1, 0))
        
        # write one x at a time
        for x in range(X):
            swizzled_dataset[w, x, :, :] = swizzled_data_w[x, :, :]

if __name__ == '__main__':
    filename = sys.argv[0]
    swizzle(filename)
