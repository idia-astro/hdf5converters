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
    def __init__(self, data, axes, swizzled=False):
        if swizzled:
            self.name = "DATA_" + axes
        else:
            self.name = "DATA"
        
        self.axes = axes        
        self.data = data
    
    # TODO axes object
    # TODO I don't think we actually need this
    #def axis_name(self, axis_numeric):
        #"""Convert numeric axes to named axes, relative to this dataset
            #e.g. if axes are XYZW, 1 -> Z, (2, 3) -> XY
        #"""
        #if isinstance(axis_numeric, int):
            #axis_numeric = (axis_numeric,)
        
        #return "".join(sorted(reversed(self.axes)[d] for d in axis_numeric))
        
    def axis_numeric(self, axis_name):
        """Convert named axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, Z -> 1, XY -> (2, 3)
        """
        axis_numeric = tuple(sorted(reversed(self.axes).index(l) for l in axis_name))
        
        if len(axis_numeric) == 1:
            return axis_numeric[0]
        
        return axis_numeric
    
    def write_statistics(self, hdu_group, axis_name):
        stats = hdu_group.require_group("Statistics/%s/%s" % (self.name, axis_name))
        
        axis = self.axis_numeric(axis_name)
        
        with warnings.catch_warnings():
            # nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stats.create_dataset("MEANS", np.nanmean(self.data, axis=axis))
            stats.create_dataset("MIN_VALS", np.nanmin(self.data, axis=axis))
            stats.create_dataset("MAX_VALS", np.nanmax(self.data, axis=axis))
            
        stats.create_dataset("NAN_COUNTS", np.count_nonzero(np.isnan(self.data), axis=axis))
        
    #def write_histogram(self, path, data, axis):
        ## TODO make this generic; able to work with any data and axis
        ## TODO CURRENTLY UNTESTED for 4D
        #dims = data.shape
        #W = dims[0]
        #Z = dims[1]
        #N = int(np.sqrt(dims[2] * dims[3]))
        
        #bins = np.zeros([W, Z, N])
        #edges = np.zeros([2, W, Z])
        
        #for j in range(W):
            #for i in range(Z):
                #data_slice = data[j, i, :, :]
                #data_notnan = data_slice[~np.isnan(data_slice)]
                #if data_notnan.size:
                    #b, e = np.histogram(data_notnan, N)
                    #bins[j, i, :] = b
                    #edges[:, j, i] = e[:2]
                #else:
                    #bins[j, i, :] = np.nan
                    #edges[:, j, i] = np.nan
        
        #widths = edges[1] - edges[0]
        #first_centers = (edges[0] + edges[1]) / 2
        
        #path = path + "/Histograms"
        
        #self.write_data(path, "BINS", bins)
        #self.write_data(path, "BIN_WIDTHS", widths)
        #self.write_data(path, "FIRST_CENTERS", first_centers)
    
    #def write_percentiles(self, path, data, axis):
        #percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])
        
        #with warnings.catch_warnings():
            ## nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            #warnings.simplefilter("ignore", category=RuntimeWarning)
            #percentile_values = np.nanpercentile(data, percentiles, axis=axis).transpose() # TODO test with other axes; see if transposing is always the right thing to do
            
        #path = path + "/Percentiles"
            
        #self.write_data(path, "PERCENTILES", percentiles)
        #self.write_data(path, "VALUES", percentile_values)
    
    def write(self, hdu_group, statistics_axes, chunks):        
        # write this dataset
        hdu_group.create_dataset(self.name, data=self.data, chunks=chunks)
        
        # write statistics
        for axis_name in statistics_axes:
            self.write_statistics(hdu_group, axis_name)
            #self.write_histogram(hdu_group, axis_name)
            #self.write_percentiles(hdu_group, axis_name)


class HDUGroup:
    # FITS header attributes to keep (exact names)
    FITS_KEEP = ('BUNIT', 'DATE-OBS', 'EQUINOX', 'INSTR', 'OBSDEC', 'OBSERVER', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z', 'OBSRA', 'RADESYS', 'TELE', 'TIMESYS')
    
    # FITS header attributes to keep (regular expression matches)
    FITS_KEEP_RE = (
        re.compile('^(CDELT|CROTA|CRPIX|CRVAL|CTYPE|CUNIT)\d+'),
    )
    
    def __init__(self, name, fits_hdu):
        self.name = name
        self.data = fits_hdu.data
        self.header = fits_hdu.header
        
    @staticmethod
    def convert(val):
        if isinstance(val, str):
            return np.string_(val)
        return val

    def copy_attrs(self, hdu_group, keys):
        for key in keys:
            if key in self.header:
                hdu_group.attrs.create(key, Converter.convert(self.header[key]))
        
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
            
        # TODO get the actual axes out of the header
        main_dataset = DataSet(self.data, "XYZW")
        
        # TODO get statistics axes from args
        main_dataset.write(hdu_group, statistics_axes=statistics_axes, chunks=chunks)
        
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
        primary = HDUGroup("Primary", self.fits[0])
        primary.write(self.hdf5, statistics_axes=("XY", "Z", "XYZ"), chunks=tuple(args.chunks) if args.chunks else None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input filename')
    # TODO how do we change this when we keep the 4D dataset?
    parser.add_argument('--chunks', nargs=3, type=int, help='Chunks to use, order: Z Y X')
    args = parser.parse_args()
    
    baseFileName, _ = os.path.splitext(args.filename)

    if args.chunks:
        outputFileName = baseFileName + "_chunked_{}_{}_{}.hdf5".format(args.chunks[0], args.chunks[1], args.chunks[2])
    else:
        outputFileName = baseFileName + ".hdf5"
        
    with Converter(args.filename, outputFileName) as converter:
        converter.convert(args)
    
