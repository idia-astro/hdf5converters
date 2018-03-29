#!/usr/bin/env python3

import h5py
import numpy as np
from mpi4py import MPI
import sys

rank = MPI.COMM_WORLD.rank
nprocs = MPI.COMM_WORLD.size

# copy of serial implementation
def serial(data, dataset):
    if rank == 0:
        dataset[:,:,:,:] = np.transpose(data, (0, 3, 2, 1))

# data -> dataset directly - SLOWEST
def one(data, dataset):
    W, Z, Y, X = data.shape
    w = rank
    
    data_w = data[w].copy() # read into memory first
        
    for z in range(Z):
        for y in range(Y):
            dataset[w, :, y, z] = data_w[z, y, :]

# via transposed view (this slice only)
def two(data, dataset):
    w = rank
    dataset[w] = np.transpose(data[w], (2, 1, 0))

# transpose in memory first via copy -- FASTEST?
def three(data, dataset):
    w = rank
    dataset[w] = np.transpose(data[w], (2, 1, 0)).copy()
    
# collective -- slightly slower
def three_collective(data, dataset):
    w = rank
    with dataset.collective:
        dataset[w] = np.transpose(data[w], (2, 1, 0)).copy()

# transpose in memory first by hand
def four(data, dataset):
    W, Z, Y, X = data.shape
    w = rank
    
    swizzled_data_w = np.zeros((X, Y, Z), dtype=data.dtype)
    
    for z in range(Z):
        for y in range(Y):
            swizzled_data_w[:, y, z] = data[w, z, y, :]
    
    dataset[w] = swizzled_data_w
    
# read into memory *and* write transpose to memory -- doesn't seem to help; slightly slower

def five(data, dataset):
    w = rank
    dataset[w] = np.transpose(data[w].copy(), (2, 1, 0)).copy()

FUNCS = {
    "serial": serial,
    "one": one,
    "two": two,
    "three": three,
    "three_collective": three_collective,
    "four": four,
    "five": five,
}


def swizzle(filename, funcname): # for this prototype assume XYZW -> ZYXW
    with h5py.File(filename, 'r+', driver='mpio', comm=MPI.COMM_WORLD) as f:
        data = f["0"]["DATA"]
        W, Z, Y, X = data.shape
        
        swizzled_group = f["0"].require_group("SwizzledData")
        swizzled_dataset = swizzled_group.create_dataset("ZYXW", (W, X, Y, Z), dtype=data.dtype)
        
        func = FUNCS[funcname]
        func(data, swizzled_dataset)


if __name__ == '__main__':
    filename, funcname = sys.argv[1:]
    swizzle(filename, funcname)
