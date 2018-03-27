#!/usr/bin/env python3

import h5py
import numpy as np
from mpi4py import MPI
import argparse
import itertools

rank = MPI.COMM_WORLD.rank
nprocs = MPI.COMM_WORLD.size

def swizzle(args): # for this prototype assume XYZW -> ZYXW
    with h5py.File(args.filename, 'r+', driver='mpio', comm=MPI.COMM_WORLD) as f:
        data = f["0"]["DATA"]
        W, Z, Y, X = data.shape
        
        swizzled_group = f["0"].require_group("SwizzledData")
        swizzled_data = swizzled_group.create_dataset("ZYXW", (W, X, Y, Z), dtype=data.dtype)
        
        w = rank
        
        for z in range(Z):
            for y in range(Y):
                swizzled_data[w, :, y, z] = data[w, z, y, :]



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tool for swizzling datasets using MPI.")
    
    parser.add_argument('filename', help='Input filename')
        
    args = parser.parse_args()
    
    swizzle(args)
