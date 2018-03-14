# hdf5converters
Python scripts to convert FITS images to HDF5 images using the [custom IDIA schema](https://github.com/idia-astro/hdf5converters/wiki/HDF5-Image-Schema).

## `fits2hdf5.py`
`fits2hdf5.py` converts a single FITS file to HDF5 using the IDIA schema. The original data and a selection of header attributes are copied, and configurable statistics and swizzled datasets may also be written. Currently only the primary HDU from the original file is processed.
