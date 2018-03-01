#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import warnings
import sys
import argparse
import os
import re

class Converter:
    # FITS header attributes to keep (exact names)
    FITS_KEEP = ('BUNIT', 'DATE-OBS', 'EQUINOX', 'INSTR', 'OBSDEC', 'OBSERVER', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z', 'OBSRA', 'RADESYS', 'TELE', 'TIMESYS')
    
    # FITS header attributes to keep (regular expression matches)
    FITS_KEEP_RE = (
        re.compile('^(CDELT|CROTA|CRPIX|CRVAL|CTYPE|CUNIT)\d+'),
    )
    
    def __init__(self, fitsname, hdf5name):
        self.fitsname = fitsname
        self.hdf5name = hdf5name
        
    def __enter__(self):
        self.fits = fits.open(self.fitsname)
        self.hdf5 = h5py.File(self.hdf5name, "w")
        
        data = self.fits[0].data
        dims = data.shape
        
        if len(dims) == 4:
            dims = dims[1:]
            data = data[0, :, :, :]
            
        elif len(dims) == 2:
            dims = (1, *dims)
            data = data[None, :]
            
        elif len(dims) == 3:
            print('3D data found.')
        
        else:
            raise ValueError('Could not coerce data to 3D cube.')
        
        print('File dims: {}'.format(dims))
        print('Chunk dims: {}'.format(args.chunks))
        
        self.data = data
        self.dims = dims
        self.header = self.fits[0].header
            
        return self
        
    def __exit__(self, e_type, e_value, e_traceback):
        del self.data
        del self.header
        self.fits.close()
        self.hdf5.close()

    @staticmethod
    def convert(val):
        if isinstance(val, str):
            return np.string_(val)
        return val

    def copy_attrs(self, path, keys):
        group = self.hdf5.require_group(path)
        
        for key in keys:
            if key in self.header:
                group.attrs.create(key, Converter.convert(self.header[key]))
    
    def write_data(self, path, name, data, chunks=None):
        group = self.hdf5.require_group(path)
        group.require_dataset(name, data=data, chunks=chunks)
        
    def write_attr(self, path, name, attr):
        group = self.hdf5.require_group(path)
        group.attrs.create(name, attr)
    
    def write_average_data(self, path, name, data, axis):
        average_data = np.nanmean(data, axis=axis)
        self.write_data(path, name, average_data)
        return average_data
    
    def write_statistics(self, path, data, axis):
        with warnings.catch_warnings():
            # nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            means = np.nanmean(data, axis=axis)
            min_vals = np.nanmin(data, axis=axis)
            max_vals = np.nanmax(data, axis=axis)
            
        nan_counts = np.count_nonzero(np.isnan(data), axis=axis)
        
        self.write_data(path, "MEANS", means)
        self.write_data(path, "MIN_VALS", min_vals)
        self.write_data(path, "MAX_VALS", max_vals)
        self.write_data(path, "NAN_COUNTS", nan_counts)
            
    def write_histogram(self, path, data, axis):
        # TODO make this generic; able to work with any data and axis
        dims = data.shape
        Z = dims[0]
        N = int(np.sqrt(dims[1] * dims[2]))
        
        bins = np.zeros([Z, N])
        edges = np.zeros([2, Z])
        for i in range(Z):
            data_slice = data[i, :, :]
            data_notnan = data_slice[~np.isnan(data_slice)]
            if data_notnan.size:
                b, e = np.histogram(data_notnan, N)
                bins[i, :] = b
                edges[:,i] = e[:2]
            else:
                bins[i, :] = np.nan
                edges[:,i] = np.nan
        
        widths = edges[1] - edges[0]
        first_centers = (edges[0] + edges[1]) / 2
        
        path = path + "/Histograms"
        
        self.write_data(path, "BINS", bins)
        self.write_data(path, "BIN_WIDTHS", widths)
        self.write_data(path, "FIRST_CENTERS", first_centers)
    
    def write_percentiles(self, path, data, axis):
        percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])
        
        with warnings.catch_warnings():
            # nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            percentile_values = np.nanpercentile(data, percentiles, axis=axis).transpose() # TODO test with other axes; see if transposing is always the right thing to do
            
        path = path + "/Percentiles"
            
        self.write_data(path, "PERCENTILES", percentiles)
        self.write_data(path, "VALUES", percentile_values)
        
    def convert(self, args):
        # Copy attributes from header
        attrs_to_copy = set(FITS_KEEP)
        for regex in FITS_KEEP_RE:
            attrs_to_copy |= {k for k in header if regex.search(k)}
        
        self.copy_attrs("Primary", attrs_to_copy)

        if 'HISTORY' in self.header:
            self.write_data('Primary', 'HISTORY', [np.string_(val) for val in header['HISTORY']])
            
        if 'COMMENT' in self.header:
            self.write_data('Primary', 'COMMENT', [np.string_(val) for val in header['COMMENT']])
        
        self.write_data("Primary", "DATA", data, tuple(args.chunks) if args.chunks else None)
        
        # STATISTICS
        # TODO: average datasets
        # TODO: pass in averages, mipmaps, swizzles to calculate on commandline
        # TODO: calculate subgroups automatically from data and axis
        # TODO: submethods should also know to write to Statistics
        
        self.write_percentiles("Primary/Statistics/DATA/XY", self.data, (1, 2))
        self.write_histogram("Primary/Statistics/DATA/XY", self.data, (1, 2))
        self.write_statistics("Primary/Statistics/DATA/XY", self.data, (1, 2))
        
        # AVERAGE DATA
        
        # TODO again, most of the path should be generated automatically
        
        average_data_Z = self.write_average_data("Primary/AverageData/DATA", "Z", self.data, 0)
        
        # TODO TODO TODO statistics for average data
        
        # TODO this needs to reuse the existing functions
        
        #averageZStatsGroup = allStatsGroup.create_group("StatisticsAverageZ")
        
        #averageDataNaNFixed = averageData[~np.isnan(averageData)]
        #means_avg_z = np.mean(averageDataNaNFixed)
        #minVals_avg_z = np.min(averageDataNaNFixed)
        #maxVals_avg_z = np.max(averageDataNaNFixed)
        
        #nanCounts_avg_z = np.count_nonzero(np.isnan(averageData))
        
        #percentileVals_avg_z = np.percentile(averageDataNaNFixed, percentiles)
        
        #(tmpBins, tmpEdges) = np.histogram(averageDataNaNFixed, N)
        #histogramBins_avg_z = tmpBins
        #histogramBinWidths_avg_z = tmpEdges[1] - tmpEdges[0]
        #histogramFirstBinCenters_avg_z = (tmpEdges[0] + tmpEdges[1]) / 2.0
        
        
        #averageZStatsGroup.create_dataset("MEANS", [1], dtype='f4', data=means_avg_z)
        #averageZStatsGroup.create_dataset("MIN_VALS", [1], dtype='f4', data=minVals_avg_z)
        #averageZStatsGroup.create_dataset("MAX_VALS", [1], dtype='f4', data=maxVals_avg_z)
        #averageZStatsGroup.create_dataset("NAN_COUNTS", [1], dtype='i4', data=nanCounts_avg_z)

        #averageZhistGroup = averageZStatsGroup.create_group("Histograms")
        #averageZhistGroup.create_dataset("FIRST_CENTERS", [1], dtype='f4', data=histogramFirstBinCenters_avg_z)
        #averageZhistGroup.create_dataset("BIN_WIDTHS", [1], dtype='f4', data=histogramBinWidths_avg_z)
        #averageZhistGroup.create_dataset("BINS", [1, N], dtype='i4', data=histogramBins_avg_z)

        #averageZpercentileGroup = averageZStatsGroup.create_group("Percentiles")
        #averageZpercentileGroup.create_dataset("PERCENTILES", [len(percentiles)], dtype='f4', data=percentiles)
        #averageZpercentileGroup.create_dataset("VALUES", [1, len(percentiles)], dtype='f4', data=percentileVals_avg_z)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input filename')
    parser.add_argument('--chunks', nargs=3, type=int, help='Chunks to use, order: Z Y X')
    parser.add_argument('--stokes', help='Stokes parameter to assign', default='')
    args = parser.parse_args()
    
    if args.stokes:
        print('Using Stokes parameter {}.'.format(args.stokes))
    else:
        print('No Stokes parameter specified.')
    
    baseFileName, _ = os.path.splitext(args.filename)

    if args.chunks:
        outputFileName = baseFileName + "_chunked_{}_{}_{}.hdf5".format(args.chunks[0], args.chunks[1], args.chunks[2])
    else:
        outputFileName = baseFileName + ".hdf5"
        
    with Converter(inputFileName, outputFileName) as converter:
        converter.convert(args)
    
