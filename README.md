# hdf5converters
Python module to convert FITS images to HDF5 images using the [custom IDIA schema](https://github.com/idia-astro/hdf5converters/wiki/HDF5-Image-Schema).

<a title="User:MagicImage, Public domain, via Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File:Stop_sign_light_red.svg"><img width="64" align="left" alt="Stop sign light red" src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/9d/Stop_sign_light_red.svg/64px-Stop_sign_light_red.svg.png"></a>

This is an outdated implementation which is not compatible with recent versions of the IDIA schema. Please use [`fits2idia`](https://github.com/CARTAvis/fits2idia) instead.

## Executable scripts

`fits2hdf5` converts a single FITS file to HDF5 using the IDIA schema. The original data and a selection of header attributes are copied, and configurable statistics and swizzled datasets may also be written. Currently only the primary HDU from the original file is processed. To see a list of commandline parameters:

    fits2hdf5 --help

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

## Defaults

By default both the function and the script are configured to generate a core set of additional datasets which is expected by the IDIA HDF5 viewer. These parameters can be overriden with custom parameters.
