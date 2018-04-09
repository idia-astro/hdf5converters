#!/bin/bash

filename=$1; shift
nprocs=$1; shift

hdf5filename="${filename%.fits}.hdf5"

function drop_caches {
    # sync before dropping caches
    sudo sync
    # wait
    sleep 5
    # drop caches
    sudo bash -c 'echo 3 > /proc/sys/vm/drop_caches'
}

# convert without swizzling
echo "Convert FITS to HDF5"
fits2hdf5 $filename

# time serial swizzling and save a serially swizzled copy to compare for correctness
cp $hdf5filename serial_reference.hdf5
echo "Time serial swizzle"

drop_caches
/usr/bin/time -v hdf5swizzle_serial.py serial_reference.hdf5

for funcname in $@
do
    echo -e "\n-------------------------------\nSwizzle with function $funcname"
    
    # recreate file
    rm -f test.hdf5
    cp $hdf5filename test.hdf5
    
    # time parallel swizzle
    drop_caches
    /usr/bin/time -v mpiexec --mca orte_base_help_aggregate 0 -n $nprocs hdf5swizzle.py test.hdf5 $funcname
    
    # check for correctness
    if `h5diff serial_reference.hdf5 test.hdf5 0/SwizzledData/ZYXW &> /dev/null`
    then 
        echo "Dataset matches reference."
    else
        echo "!!! Alert: dataset does not match reference."
    fi
done

rm serial_reference.hdf5  test.hdf5
