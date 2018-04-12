#!/usr/bin/env python3

import h5py
import numpy as np
from mpi4py import MPI
import sys
import itertools

comm = MPI.COMM_WORLD
rank = comm.rank
nprocs = comm.size

class Swizzler:
    @staticmethod
    def w_nocopy(data, dataset):
        w = rank
        dataset[w] = np.transpose(data[w], (2, 1, 0))

    @staticmethod
    def w_copy(data, dataset):
        w = rank
        dataset[w] = np.transpose(data[w], (2, 1, 0)).copy()

    @staticmethod
    def wz_nocopy(data, dataset):
        W, Z, Y, X = data.shape
        w = rank
        
        data_w = np.empty((Z, Y, X), dtype=data.dtype)
        
        # read one z at a time
        for z in range(Z):
            data_w[z] = data[w,z,:,:].copy()
            
        # swizzle whole w by transposing
        swizzled_data_w = np.transpose(data_w, (2, 1, 0))
        
        # write one x at a time
        for x in range(X):
            dataset[w, x, :, :] = swizzled_data_w[x, :, :]

    @staticmethod
    def wz_copy(data, dataset):
        W, Z, Y, X = data.shape
        w = rank
        
        data_w = np.empty((Z, Y, X), dtype=data.dtype)
        
        # read one z at a time
        for z in range(Z):
            data_w[z] = data[w,z,:,:].copy()
            
        # swizzle whole w by transposing and copying
        swizzled_data_w = np.transpose(data_w, (2, 1, 0)).copy()
        
        # write one x at a time
        for x in range(X):
            dataset[w, x, :, :] = swizzled_data_w[x, :, :]


def swizzle(filename, funcname): # for this prototype assume XYZW -> ZYXW
    with h5py.File(filename, 'r+', driver='mpio', comm=MPI.COMM_WORLD) as f:
        data = f["0"]["DATA"]
        W, Z, Y, X = data.shape
        
        swizzled_group = f["0"].require_group("SwizzledData")
        swizzled_dataset = swizzled_group.create_dataset("ZYXW", (W, X, Y, Z), dtype=data.dtype)
        
        func = getattr(Swizzler, funcname)
        func(data, swizzled_dataset)


if __name__ == '__main__':
    filename, funcname = sys.argv[1:]
    swizzle(filename, funcname)
