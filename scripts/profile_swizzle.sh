#!/bin/bash

filename=$1
nprocs=$2

hdf5filename="${filename%.fits}.hdf5"

function drop_caches {
    # sync before dropping caches
    sudo sync
    # wait
    sleep 5
    # drop caches
    sudo bash -c 'echo 3 > /proc/sys/vm/drop_caches'
}

drop_caches

# convert with swizzled data and time
echo "Serial swizzle (includes write of original dataset)"
/usr/bin/time -v fits2hdf5 --swizzles ZYXW -- $filename

# save a serially swizzled copy to compare for correctness
mv $hdf5filename serial_reference.hdf5

# convert without swizzled data
fits2hdf5 $filename

for funcname in "serial" "one" "two" "three" "three_collective" "four" "five"
do
    echo -e "\n-------------------------------\nSwizzle with function $funcname"
    
    # recreate file
    rm -f test.hdf5
    cp $hdf5filename test.hdf5
    
    # time parallel swizzle
    drop_caches
    /usr/bin/time -v mpiexec -n $nprocs hdf5swizzle.py test.hdf5 $funcname
    
    # check for correctness
    if `h5diff serial_reference.hdf5 test.hdf5 0/SwizzledData/ZYXW &> /dev/null`
    then 
        echo "Dataset matches reference."
    else
        echo "!!! Alert: dataset does not match reference."
    fi
done

rm serial_reference.hdf5  test.hdf5
