#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import warnings
import argparse
import os
import re

# helper class storing state for an original or swizzled dataset
class DataSet:
    def __init__(self, hdu_group, data, axes, swizzled=False):
        self.hdu_group = hdu_group
        
        if swizzled:
            self.name = "DATA_" + axes
        else:
            self.name = "DATA"
        
        self.axes = axes
        self._reverse_axes = list(reversed(axes))
        self.data = data
    
    # TODO axes object
    # TODO I don't think we actually need this
    #def axis_name(self, axis_numeric):
        #"""Convert numeric axes to named axes, relative to this dataset
            #e.g. if axes are XYZW, 1 -> Z, (2, 3) -> XY
        #"""
        #if isinstance(axis_numeric, int):
            #axis_numeric = (axis_numeric,)
        
        #return "".join(sorted(self._reverse_axes[d] for d in axis_numeric))
        
    def axis_numeric(self, axis_name):
        """Convert named axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, Z -> 1, XY -> (2, 3)
        """
        axis_numeric = tuple(sorted(self._reverse_axes.index(l) for l in axis_name))
        
        if len(axis_numeric) == 1:
            return axis_numeric[0]
        
        return axis_numeric
    
    def write_statistics(self, axis_name):
        stats = self.hdu_group.require_group("Statistics/%s/%s" % (self.name, axis_name))
        axis = self.axis_numeric(axis_name)
        
        with warnings.catch_warnings():
            # nanmean, etc. print a warning for empty slices, i.e. planes full of nans, but give the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stats.create_dataset("MEAN", data=np.nanmean(self.data, axis=axis))
            stats.create_dataset("MIN", data=np.nanmin(self.data, axis=axis))
            stats.create_dataset("MAX", data=np.nanmax(self.data, axis=axis))
            
        stats.create_dataset("NAN_COUNT", data=np.count_nonzero(np.isnan(self.data), axis=axis))
        
    def write_histogram(self, axis_name):
        stats = self.hdu_group.require_group("Statistics/%s/%s" % (self.name, axis_name))
        axis = self.axis_numeric(axis_name)
        
        dims = data.shape
        N = np.sqrt(np.multiply.reduce([dims[a] for a in axis]))
        
        # TODO now use a recursive function (?)
                
        #dims = data.shape
        #W = dims[0]
        #Z = dims[1]
        #N = int(np.sqrt(dims[2] * dims[3]))
        
        #bins = np.zeros([W, Z, N])
        
        #for j in range(W):
            #for i in range(Z):
                #data_slice = data[j, i, :, :]
                #data_notnan = data_slice[~np.isnan(data_slice)]
                #if data_notnan.size:
                    #b, _ = np.histogram(data_notnan, N)
                    #bins[j, i, :] = b
                #else:
                    #bins[j, i, :] = np.nan


        stats.create_dataset("HISTOGRAM", bins)
    
    #def write_percentiles(self, path, data, axis):
        percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])
        self.hdu_group.create_dataset("PERCENTILE_RANKS", percentiles)
        
        stats = self.hdu_group.require_group("Statistics/%s/%s" % (self.name, axis_name))
        
        # TODO rewrite this
        
        #with warnings.catch_warnings():
            ## nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            #warnings.simplefilter("ignore", category=RuntimeWarning)
            #percentile_values = np.nanpercentile(data, percentiles, axis=axis).transpose() # TODO test with other axes; see if transposing is always the right thing to do
        
        stats.create_dataset("PERCENTILES", percentile_values)
    
    def write(self, statistics_axes, chunks):        
        # write this dataset
        self.hdu_group.create_dataset(self.name, data=self.data, chunks=tuple(chunks) if chunks else None)
        
        # write statistics
        for axis_name in statistics_axes:
            self.write_statistics(axis_name)
            #self.write_histogram(axis_name)
            #self.write_percentiles(axis_name)


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
        
    def write(self, hdf5file, statistics_axes=None, chunks=None):
        hdu_group = hdf5file.require_group(self.name)
        
        # Copy attributes from header
        attrs_to_copy = set(self.FITS_KEEP)
        for regex in self.FITS_KEEP_RE:
            attrs_to_copy |= {k for k in self.header if regex.search(k)}
        
        self.copy_attrs(hdu_group, attrs_to_copy)

        if 'HISTORY' in self.header:
            hdu_group.create_dataset('HISTORY', [np.string_(val) for val in self.header['HISTORY']])
            
        if 'COMMENT' in self.header:
            hdu_group.create_dataset('COMMENT', [np.string_(val) for val in self.header['COMMENT']])
            
        num_axes = self.header["NAXIS"]
        
        # TODO: this is good enough for now, but we should look for a more robust way to detect the order
        axes = "XYZW"[:num_axes]
        
        main_dataset = DataSet(hdu_group, self.data, axes)
        
        main_dataset.write(statistics_axes=statistics_axes, chunks=chunks)
        
        # TODO: swizzles


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
        # TODO: iterate over all the HDUs -- how do we pick the names?
        # TODO: some HDUs are just tables.
        primary = HDUGroup("Primary", self.fits[0])
        primary.write(self.hdf5, statistics_axes=args.statistics_axes, chunks=args.chunks)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input filename')
    # TODO how do we change this when we keep the 4D dataset?
    parser.add_argument('--chunks', nargs=3, type=int, help='Chunks to use, order: Z Y X')
    parser.add_argument('--statistics-axes', nargs="+", help='Axes along which statistics should be calculated, e.g. XY, Z, XYZ')
    # TODO if we want to split out stokes, we should pass in a stokes parameter
    # TODO optional parameter to override axes?
    args = parser.parse_args()
    
    baseFileName, _ = os.path.splitext(args.filename)

    if args.chunks:
        outputFileName = baseFileName + "_chunked_{}_{}_{}.hdf5".format(args.chunks[0], args.chunks[1], args.chunks[2])
    else:
        outputFileName = baseFileName + ".hdf5"
        
    with Converter(args.filename, outputFileName) as converter:
        converter.convert(args)
    
