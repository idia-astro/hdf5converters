#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import warnings
import sys
import argparse
import os
import re

# TODO: helper class storing state for an original or swizzled dataset
# Parent class should not be called directly
class Data:
    def __init__(self, axes, average_axes, data, chunks):
        self.path = None # not implemented
        self.stats_path = None # not implemented
        
        self.axes = axes
        self.average_axes = average_axes
        
        self.data = data
        
    def axis_name(axis_numeric):
        """ Convert numeric axes to named axes, relative to this dataset
            e.g. if axes are XYZW, 1 -> Z, (2, 3) -> XY
        """
        if isinstance(axis_numeric, int):
            axis_numeric = (axis_numeric,)
        
        return "".join(sorted(reversed(self.axes)[d] for d in axis_numeric))
        
    def axis_numeric(axis_name):
        """ Convert named axes to numeric axes, relative to this dataset
            e.g. if axes are XYZW, Z -> 1, XY -> (2, 3)
        """
        axis_numeric = tuple(sorted(reversed(self.axes).index(l) for l in axis_name))
        
        if len(axis_numeric) == 1:
            return axis_numeric[0]
        
        return axis_numeric
    
    def write_statistics(self, parent, data, axis_name):
        axis = self.axis_numeric(axis_name)
        
        with warnings.catch_warnings():
            # nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            warnings.simplefilter("ignore", category=RuntimeWarning)
            parent.create_dataset("MEANS", np.nanmean(data, axis=axis))
            parent.create_dataset("MIN_VALS", np.nanmin(data, axis=axis))
            parent.create_dataset("MAX_VALS", np.nanmax(data, axis=axis))
            
        parent.create_dataset("NAN_COUNTS", np.count_nonzero(np.isnan(data), axis=axis))
    
    def write(self, parent, statistics_axes=None, chunks=None):
        parent_group = self.hdf5.require_group(parent)
        
        # write this dataset
        group = parent.require_group(self.path)
        group.create_dataset(self.path, data=self.data, chunks=chunks)
        
        for axis_name in statistics_axes:
            # write statistics
            stats = parent.require_group(self.stats_path)
            
            self.write_statistics(stats, data, axis_name)
            
            # histogram
            # percentiles
            
    def write_average(self, parent, average_axis, statistics_axes=None):
        # we need to pass in the mean data
        # but then we need to collapse the axes -- TODO!
        average_axis_numeric = self.axis_numeric(average_axis)
        average_data = AverageData(self.axes, average_axes, np.nanmean(self.data, axis=average_axis_numeric))
        average_data.write(parent, statistics_axes) 

# This can be called directly from the converter
class OriginalData(Data):
    def __init__(self, data):
        super(OriginalData, self).__init__("XYZW", None, data)
        self.path = "DATA"
        self.stats_path = "Statistics/%s" % self.path

# This can be called directly from the converter
class SwizzledData(Data):
    def __init__(self, axes, data):
        super(SwizzledData, self).__init__("SwizzledData", "DATA_%s" % axes, axes, None, data)
        self.path = "DATA_%s" % axes
        self.stats_path = "Statistics/%s" % self.path
        
# This should be only created from the original or swizzled dataset object, not called separately from the constructor
class AverageData(Data):
    def __init__(self, axes, average_axes, data):
        super(AverageData, self).__init__(axes, average_axes, data)
        if axes == "XYZW": # average of original
            source = "DATA"
        else: # average of swizzled
            source = "DATA_%s" % axes
            
        self.path = "AverageData/%s/%s" % (source, average_axes)
        self.stats_path = "Statistics/%s_AVG_%s"  % (source, average_axes)

# TODO decide on 3D vs 4D

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
            print('4D data found.')
            
        elif len(dims) == 3:
            print('3D data found.')
            dims = (1, *dims)
            data = data[None, :]
            
        elif len(dims) == 2:
            print('2D data found.')
            dims = (1, 1, *dims)
            data = data[None, None :]
        
        else:
            raise ValueError('Could not coerce data to 4D.')
        
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
        
    def write_attr(self, path, name, attr):
        group = self.hdf5.require_group(path)
        group.attrs.create(name, attr)
    
    #def write_average_data(self, path, name, data, axis):
        #average_data = np.nanmean(data, axis=axis)
        #self.write_data(path, name, average_data)
        #return average_data
    
    #def write_statistics(self, path, data, axis):
        #with warnings.catch_warnings():
            ## nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
            #warnings.simplefilter("ignore", category=RuntimeWarning)
            #means = np.nanmean(data, axis=axis)
            #min_vals = np.nanmin(data, axis=axis)
            #max_vals = np.nanmax(data, axis=axis)
            
        #nan_counts = np.count_nonzero(np.isnan(data), axis=axis)
        
        #self.write_data(path, "MEANS", means)
        #self.write_data(path, "MIN_VALS", min_vals)
        #self.write_data(path, "MAX_VALS", max_vals)
        #self.write_data(path, "NAN_COUNTS", nan_counts)
            
    def write_histogram(self, path, data, axis):
        # TODO make this generic; able to work with any data and axis
        # TODO CURRENTLY UNTESTED for 4D
        dims = data.shape
        W = dims[0]
        Z = dims[1]
        N = int(np.sqrt(dims[2] * dims[3]))
        
        bins = np.zeros([W, Z, N])
        edges = np.zeros([2, W, Z])
        
        for j in range(W):
            for i in range(Z):
                data_slice = data[j, i, :, :]
                data_notnan = data_slice[~np.isnan(data_slice)]
                if data_notnan.size:
                    b, e = np.histogram(data_notnan, N)
                    bins[j, i, :] = b
                    edges[:, j, i] = e[:2]
                else:
                    bins[j, i, :] = np.nan
                    edges[:, j, i] = np.nan
        
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
        # TODO: axes should be calculated automatically
        
        self.write_percentiles("Primary/Statistics/DATA/XY", self.data, (2, 3))
        self.write_histogram("Primary/Statistics/DATA/XY", self.data, (2, 3))
        self.write_statistics("Primary/Statistics/DATA/XY", self.data, (2, 3))
        
        # AVERAGE DATA
        
        # TODO again, most of the path should be generated automatically
        
        average_data_Z = self.write_average_data("Primary/AverageData/DATA", "Z", self.data, 1)
        
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
    
