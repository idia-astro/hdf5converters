#!/usr/bin/env python3

import h5py
import numpy as np
import sys

# for this prototype assume XYZW -> ZYXW
if __name__ == '__main__':
    filename = sys.argv[1]
    with h5py.File(filename, 'r+') as f:
        data = f["0"]["DATA"]
        swizzled_group = f["0"].require_group("SwizzledData")
        swizzled_dataset = swizzled_group.create_dataset("ZYXW", data=np.transpose(data, (0, 3, 2, 1)).copy())
