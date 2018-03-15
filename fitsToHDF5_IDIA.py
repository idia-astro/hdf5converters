#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import warnings
import argparse
import os
import re
import itertools
import logging


# helper class storing state for an original or swizzled dataset
class Dataset:
    def __init__(self, hdu_group, data, axes):
        self.hdu_group = hdu_group
        self.name = "DATA"
        self.axes = axes
        self._reverse_axes = list(reversed(axes))
        self.data = data
        
    def axis_numeric(self, axis_name):
        """Convert named axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, Z -> 1, XY -> (2, 3)
        """
        return tuple(sorted(self._reverse_axes.index(l) for l in axis_name))
        
    def swizzle_axis_numeric(self, axis_name):
        """Convert named swizzle axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, ZYXW -> (0, 3, 2, 1)
        """
        return tuple(self._reverse_axes.index(l) for l in reversed(axis_name))
    
    def max_bins(self):
        """Limit the number of histogram bins to the root of the size of the largest product of two dimensions.
            This will probably be the root of X * Y.
        """
        return int(np.sqrt(max(a * b for a, b in itertools.combinations(self.data.shape, 2))))
    
    def write_statistics(self, axis_name):
        if not all(d in self.axes for d in axis_name):
            logging.warning("Could not average %s dataset along %s." % (self.axes, axis_name))
            return
        
        axis = self.axis_numeric(axis_name)
        
        data_size = np.multiply.reduce([self.data.shape[d] for d in axis])
        
        if data_size == 1:
            logging.warning("Not calculating statistics for %s dataset averaged along %s: data size of 1." % (self.axes, axis_name))
            return
        
        logging.info("Writing statistics for axis %s..." % axis_name)
        
        stats = self.hdu_group.require_group("Statistics/%s" % axis_name)
        
        with warnings.catch_warnings():
            # nanmean, etc. print a warning for empty slices, i.e. planes full of nans, but give the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stats.create_dataset("MEAN", data=np.nanmean(self.data, axis=axis))
            stats.create_dataset("MIN", data=np.nanmin(self.data, axis=axis))
            stats.create_dataset("MAX", data=np.nanmax(self.data, axis=axis))
            
        stats.create_dataset("NAN_COUNT", data=np.count_nonzero(np.isnan(self.data), axis=axis))
        
    def write_histogram(self, axis_name):
        if not all(d in self.axes for d in axis_name):
            logging.warning("Could not average %s dataset along %s." % (self.axes, axis_name))
            return
        
        axis = self.axis_numeric(axis_name)
        
        shape = self.data.shape
        ndim = self.data.ndim
        
        data_size = np.multiply.reduce([shape[d] for d in axis])
        
        if data_size == 1:
            logging.warning("Not calculating histogram for %s dataset averaged along %s: data size of 1." % (self.axes, axis_name))
            return
        
        logging.info("Writing histograms for axis %s..." % axis_name)
        
        stats = self.hdu_group.require_group("Statistics/%s" % axis_name)
        
        N = min(int(np.sqrt(data_size)), self.max_bins())
                
        not_axis = [d for d in range(ndim) if d not in axis]
        data_index = [slice(None)] * ndim
        
        bins = np.zeros([shape[d] for d in not_axis] + [N])
        bin_index = [slice(None)] * bins.ndim
        
        for iterdims in itertools.product(*(range(shape[d]) for d in not_axis)):
            for d, v in zip(not_axis, iterdims):
                data_index[d] = v
                
            for d, v in enumerate(iterdims):
                bin_index[d] = v
                    
            data_slice = self.data[tuple(data_index)]
            data_notnan = data_slice[~np.isnan(data_slice)]
            
            if data_notnan.size:
                b, _ = np.histogram(data_notnan, N)
                bins[tuple(bin_index)] = b
            else:
                bins[tuple(bin_index)] = np.nan

        stats.create_dataset("HISTOGRAM", data=bins, dtype='int64')
    
    def write_percentiles(self, axis_name):
        if not all(d in self.axes for d in axis_name):
            logging.warning("Could not average %s dataset along %s." % (self.axes, axis_name))
            return
        
        percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])
        
        axis = self.axis_numeric(axis_name)
        
        data_size = np.multiply.reduce([self.data.shape[d] for d in axis])
        
        if data_size == 1:
            logging.warning("Not calculating percentiles for %s dataset averaged along %s: data size of 1." % (self.axes, axis_name))
            return
        
        logging.info("Writing percentiles for axis %s..." % axis_name)
        
        if "PERCENTILE_RANKS" not in self.hdu_group:
            self.hdu_group.create_dataset("PERCENTILE_RANKS", data=percentiles, dtype='float32')
        
        stats = self.hdu_group.require_group("Statistics/%s" % axis_name)
        
        with warnings.catch_warnings():
            # nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            percentile_values = np.nanpercentile(self.data, percentiles, axis=axis)
        
        stats.create_dataset("PERCENTILES", data=np.transpose(percentile_values, list(range(percentile_values.ndim))[1:] + [0]), dtype='float32')
        
    def write_swizzled_dataset(self, axis_name):
        if len(axis_name) != len(self.axes) or not all(d in self.axes for d in axis_name):
            logging.warning("Could not swizzle %s dataset to %s." % (self.axes, axis_name))
            return
        
        logging.info("Writing swizzled dataset %s..." % axis_name)
        
        axis = self.swizzle_axis_numeric(axis_name)
        swizzled = self.hdu_group.require_group("SwizzledData")
        swizzled.create_dataset(axis_name, dtype='f4', data=np.transpose(self.data, axis))
    
    def write(self, args):
        # write this dataset
        # TODO TODO TODO check that the number of chunks matches the data dimensions
        logging.info("Writing main dataset...")
        self.hdu_group.create_dataset(self.name, dtype='f4', data=self.data, chunks=tuple(args.chunks) if args.chunks else None)
        
        # write statistics
        for axis_name in args.statistics:
            self.write_statistics(axis_name)
            
        for axis_name in args.histograms:
            self.write_histogram(axis_name)
            
        for axis_name in args.percentiles:
            self.write_percentiles(axis_name)
        
        # write swizzled datasets
        for axis_name in args.swizzles:
            self.write_swizzled_dataset(axis_name)
            

class HDUGroup:
    # FITS header attributes to keep (exact names)
    FITS_KEEP = ('BUNIT', 'DATE-OBS', 'EQUINOX', 'INSTR', 'OBSDEC', 'OBSERVER', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z', 'OBSRA', 'RADESYS', 'TELE', 'TIMESYS')
    
    # FITS header attributes to keep (regular expression matches)
    FITS_KEEP_RE = (
        re.compile('^(CDELT|CROTA|CRPIX|CRVAL|CTYPE|CUNIT)\d+'),
        re.compile('^NAXIS\d*')
    )
    
    def __init__(self, hdf5file, name, fits_hdu):
        self.hdf5file = hdf5file
        self.name = name
        self.data = fits_hdu.data
        self.header = fits_hdu.header

    def copy_attrs(self, hdu_group, keys):
        def _convert(val):
            if isinstance(val, str):
                return np.string_(val)
            return val
        
        for key in keys:
            if key in self.header:
                hdu_group.attrs.create(key, _convert(self.header[key]))
        
    def write(self, args):
        hdu_group = self.hdf5file.require_group(self.name)
        
        # Copy attributes from header
        attrs_to_copy = set(self.FITS_KEEP)
        for regex in self.FITS_KEEP_RE:
            attrs_to_copy |= {k for k in self.header if regex.search(k)}
        
        self.copy_attrs(hdu_group, attrs_to_copy)

        if 'HISTORY' in self.header:
            hdu_group.create_dataset('HISTORY', data=[np.string_(val) for val in self.header['HISTORY']])
            
        if 'COMMENT' in self.header:
            hdu_group.create_dataset('COMMENT', data=[np.string_(val) for val in self.header['COMMENT']])
            
        num_axes = self.header["NAXIS"]
        
        # TODO: this is good enough for now, but we should look for a more robust way to detect the order
        axes = "XYZW"[:num_axes]
        
        main_dataset = Dataset(hdu_group, self.data, axes)
        main_dataset.write(args)


class Converter:
    def __init__(self, fitsname, hdf5name):
        self.fitsname = fitsname
        self.hdf5name = hdf5name
        
    def __enter__(self):
        self.fits = fits.open(self.fitsname)
        self.hdf5 = h5py.File(self.hdf5name, "w")
            
        return self
        
    def __exit__(self, e_type, e_value, e_traceback):
        self.fits.close()
        self.hdf5.close()

    def convert(self, args):
        # TODO: some HDUs are just tables.
        # TODO: what are our actual use cases for multiple HDUs? Do we want to apply the same arguments to all of them?
        primary = HDUGroup(self.hdf5, "0", self.fits[0])
        primary.write(args)
        
        
def convert(args):
    if args.quiet:
        logging.basicConfig(level=logging.CRITICAL)
    else:
        logging.basicConfig(level=logging.DEBUG)
    
    filedir, filename = os.path.split(args.filename)
    basefilename, _ = os.path.splitext(filename)

    output_dir = args.output_dir if args.output_dir else filedir
    output_filepath = os.path.join(output_dir, basefilename + ".hdf5")
        
    with Converter(args.filename, output_filepath) as converter:
        converter.convert(args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input filename')
    parser.add_argument('--chunks', nargs="+", type=int, help='Chunks to use, order: Z Y X')
    parser.add_argument('--statistics', nargs="+", help='Axes along which statistics should be calculated, e.g. XY, Z, XYZ', default=tuple())
    parser.add_argument('--histograms', nargs="+", help='Axes along which histograms should be calculated, e.g. XY, Z, XYZ', default=tuple())
    parser.add_argument('--percentiles', nargs="+", help='Axes along which percentiles should be calculated, e.g. XY, Z, XYZ', default=tuple())
    parser.add_argument('--swizzles', nargs="+", help='Axes for which swizzled datasets should be written, e.g. ZYXW', default=tuple())
    parser.add_argument('--output-dir', help="Output directory. By default, the directory of the original file will be used.")
    parser.add_argument('--quiet', action='store_true', help="Suppress all print output.")
    # TODO if we want to split out stokes, we should pass in a stokes parameter
    args = parser.parse_args()
    
    convert(args)
