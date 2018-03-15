# hdf5converters
Python module to convert FITS images to HDF5 images using the [custom IDIA schema](https://github.com/idia-astro/hdf5converters/wiki/HDF5-Image-Schema).

## Executable scripts

`fits2hdf5` converts a single FITS file to HDF5 using the IDIA schema. The original data and a selection of header attributes are copied, and configurable statistics and swizzled datasets may also be written. Currently only the primary HDU from the original file is processed.

## Module

To use this converter from a Python script or notebook, use the `hdf5_convert` helper function, which automatically passes certain default parameters to the converter. The only required argument to this function is the filename:

    from idiahdf5 import hdf5_convert
    hdf5_convert("yourfilenamehere.fits")

The default parameters can be overridden with keyword arguments passed after the filename. For example, you can change the output directory (by default the files will be written to the same directory as the original files):

    from idiahdf5 import hdf5_convert
    hdf5_convert("yourfilenamehere.fits", output_dir="/path/to/some/dir/here")

You can also bypass the helper function and call the underlying converter function directly, but you are then responsible for constructing the full argument object expected by the function (which should be an `argparse.Namespace` object or another object with a compatible interface):

    from idiahdf5.fits2hdf5 import convert
    # construct the args object however you like
    convert(args)
