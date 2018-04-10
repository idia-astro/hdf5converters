#!/usr/bin/env python3

import h5py
import numpy as np
import sys

class Swizzler:
    # serial implementation -- no copy
    @staticmethod
    def no_copy(data, dataset):
        dataset[:,:,:,:] = np.transpose(data, (0, 3, 2, 1))
    
    # serial implementation -- with copy
    @staticmethod
    def with_copy(data, dataset):
        dataset[:,:,:,:] = np.transpose(data, (0, 3, 2, 1)).copy()


def swizzle(filename, funcname): # for this prototype assume XYZW -> ZYXW
    with h5py.File(filename, 'r+') as f:
        data = f["0"]["DATA"]
        W, Z, Y, X = data.shape
        
        swizzled_group = f["0"].require_group("SwizzledData")
        swizzled_dataset = swizzled_group.create_dataset("ZYXW", (W, X, Y, Z), dtype=data.dtype)
        
        func = getattr(Swizzler, funcname)
        func(data, swizzled_dataset)


if __name__ == '__main__':
    filename, funcname = sys.argv[1:]
    swizzle(filename, funcname)
