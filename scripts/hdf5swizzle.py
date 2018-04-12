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
    # copy of serial implementation
    @staticmethod
    def serial(data, dataset):
        if rank == 0:
            dataset[:,:,:,:] = np.transpose(data, (0, 3, 2, 1))

    # data -> dataset directly - SLOWEST
    @staticmethod
    def one(data, dataset):
        W, Z, Y, X = data.shape
        w = rank
        
        data_w = data[w].copy() # read into memory first
            
        for z in range(Z):
            for y in range(Y):
                dataset[w, :, y, z] = data_w[z, y, :]

    # via transposed view (this slice only)
    @staticmethod
    def two(data, dataset):
        w = rank
        dataset[w] = np.transpose(data[w], (2, 1, 0))

    # transpose in memory first via copy -- FASTEST?
    @staticmethod
    def three(data, dataset):
        w = rank
        dataset[w] = np.transpose(data[w], (2, 1, 0)).copy()

    # transpose in memory first by hand
    @staticmethod
    def four(data, dataset):
        W, Z, Y, X = data.shape
        w = rank
        
        swizzled_data_w = np.empty((X, Y, Z), dtype=data.dtype)
        
        for z in range(Z):
            for y in range(Y):
                swizzled_data_w[:, y, z] = data[w, z, y, :]
        
        dataset[w] = swizzled_data_w
        
    # read into memory *and* write transpose to memory -- doesn't seem to help; slightly slower
    @staticmethod
    def five(data, dataset):
        w = rank
        dataset[w] = np.transpose(data[w].copy(), (2, 1, 0)).copy()
        
    # trying to get around the 2G limitation
    @staticmethod
    def six(data, dataset):
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
            
    # read and write individual channels; each process reads contiguous pieces
    @staticmethod
    def seven(data, dataset):
        W, Z, Y, X = data.shape
        
        num_channels = W * Z
        per_proc = int(np.ceil(num_channels / nprocs))
        
        range_min = rank * per_proc
        range_max = min((rank + 1) * per_proc, W * Z)
                
        data_wz = np.empty((range_max - range_min, Y, X), dtype=data.dtype)
        
        W_Z = list(itertools.product(range(W), range(Z)))
        
        for i in range(len(data_wz)):
            w, z = W_Z[range_min + i]
            data_wz[i, :, :] = data[w, z, :, :]
        
        swizzled_data_wz = np.transpose(data_wz, (2, 1, 0)).copy()
        
        for i in range(len(data_wz)):
            w, z = W_Z[range_min + i]
            dataset[w, :, :, z] = swizzled_data_wz[:, :, i]


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
