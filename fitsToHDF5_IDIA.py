#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import warnings
import sys
from typing import Sequence, Callable
import argparse
import os
import re

# FITS header attributes to keep (exact names)
FITS_KEEP = ('BUNIT', 'DATE-OBS', 'EQUINOX', 'INSTR', 'OBSDEC', 'OBSERVER', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z', 'OBSRA', 'RADESYS', 'TELE', 'TIMESYS')
# FITS header attributes to keep (regular expression matches)
FITS_KEEP_RE = (
    re.compile('^(CDELT|CROTA|CRPIX|CRVAL|CTYPE|CUNIT)\d+'),
)

def convert(val):
    if isinstance(val, str):
        return np.string_(val)
    return val


def get_attr_copier(header: fits.Header):
    def copy_attrs(obj, keys: Sequence[str]):
        for key in keys:
            if key in header:
                obj.attrs.create(key, convert(header[key]))
    return copy_attrs


def get_or_create_group(parent: h5py.Group, name: str):
    if name in parent:
        return parent[name]
    return parent.create_group(name)


def write_core(header: fits.Header, data: np.ndarray, outputHDF5: h5py.File, args: argparse.Namespace):
    copy_attrs = get_attr_copier(header)
    
    attrs_to_copy = set(FITS_KEEP)
    for regex in FITS_KEEP_RE:
        attrs_to_copy |= {k for k in header if regex.search(k)}
        
    # TODO: pass this in as a parameter
    hduGroup = get_or_create_group(outputHDF5, "Primary")
    
    copy_attrs(hduGroup, attrs_to_copy)

    # TODO: turn into datasets
    if 'HISTORY' in header:
        hduGroup.attrs.create('HISTORY', [np.string_(val) for val in header['HISTORY']])
        
    if 'COMMENT' in header:
        hduGroup.attrs.create('COMMENT', [np.string_(val) for val in header['COMMENT']])
    
    if args.chunks:
        dataSet = hduGroup.create_dataset("DATA", dims, dtype='f4', data=data, chunks=tuple(args.chunks))
    else:
        dataSet = hduGroup.create_dataset("DATA", dims, dtype='f4', data=data)

# TODO more OOP and factor out common stats writing code

def get_statistics(data):
    dims = data.shape
    Z = dims[0]
    N = int(np.sqrt(dims[1] * dims[2]))
    
    percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])

    # CALCULATE

    #means = np.zeros(Z)
    #minVals = np.zeros(Z)
    #maxVals = np.zeros(Z)
    nanCounts = np.zeros(Z)

    #averageData = np.zeros(dims[1:])
    averageCount = np.zeros(dims[1:])

    #percentileVals = np.zeros([Z, len(percentiles)])
    
    histogramBins = np.zeros([Z, N])
    histogramFirstBinCenters = np.zeros(Z)
    histogramBinWidths = np.zeros(Z)

    for i in range(Z):
        tmpData = data[i, :, :]
        nanArray = np.isnan(tmpData)
        nanCounts[i] = np.count_nonzero(nanArray)
        tmpDataNanFixed = tmpData[~nanArray]
        if nanCounts[i] == tmpData.shape[0] * tmpData.shape[1]:
            #means[i] = np.NaN
            #minVals[i] = np.NaN
            #maxVals[i] = np.NaN
            #percentileVals[i, :] = np.NaN
            histogramBins[i, :] = np.NaN
            histogramBinWidths[i] = np.NaN
            histogramFirstBinCenters[i] = np.NaN
        else:
            #means[i] = np.mean(tmpDataNanFixed)
            #minVals[i] = np.min(tmpDataNanFixed)
            #maxVals[i] = np.max(tmpDataNanFixed)
            #percentileVals[i, :] = np.percentile(tmpDataNanFixed, percentiles)
            (tmpBins, tmpEdges) = np.histogram(tmpDataNanFixed, N)
            histogramBins[i, :] = tmpBins
            histogramBinWidths[i] = tmpEdges[1] - tmpEdges[0]
            histogramFirstBinCenters[i] = (tmpEdges[0] + tmpEdges[1]) / 2.0

            #averageCount += (~nanArray).astype(int)
            #averageData += np.nan_to_num(tmpData)
            
    #averageData /= np.fmax(averageCount, 1)
    #averageData[averageCount < 1] = np.NaN

    # NEW VECTORISED CODE
    
    print(histogramBins)
    print(histogramBinWidths)
    print(histogramFirstBinCenters)
        
    # TODO no vectorised histogram function; can we do something clever?
        
    with warnings.catch_warnings():
        # nanmean prints a warning for empty slices, i.e. planes full of nans, but gives the correct answer (nan).
        warnings.simplefilter("ignore", category=RuntimeWarning)
        means = np.nanmean(data, axis=(1, 2))
        minVals = np.nanmin(data, axis=(1, 2))
        maxVals = np.nanmax(data, axis=(1, 2))
        percentileVals = np.nanpercentile(data, percentiles, axis=(1, 2)).transpose() # TODO test with other axes; see if transposing is always the right thing to do
        
    nanCounts = np.count_nonzero(np.isnan(data), axis=(1, 2))
    averageData = np.nanmean(data, axis=0)
                    
    return means, minVals, maxVals, nanCounts, averageData, averageCount, percentileVals, histogramBins, histogramFirstBinCenters, histogramBinWidths
    

def write_statistics(header: fits.Header, data: np.ndarray, outputHDF5: h5py.File, args: argparse.Namespace):
    dims = data.shape
    Z = dims[0]
    N = int(np.sqrt(dims[1] * dims[2]))
    
    percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])

    # CALCULATE
            
    means, minVals, maxVals, nanCounts, averageData, averageCount, percentileVals, histogramBins, histogramFirstBinCenters, histogramBinWidths = get_statistics(data)
    
    # WRITE STATISTICS FOR MAIN DATASET
    
    # TODO: pass this in as a parameter
    hduGroup = get_or_create_group(outputHDF5, "Primary")
        
    allStatsGroup = hduGroup.create_group("Statistics")
    statsGroup = allStatsGroup.create_group("StatisticsXY")
    
    statsGroup.create_dataset("MEANS", [Z], dtype='f4', data=means)
    statsGroup.create_dataset("MIN_VALS", [Z], dtype='f4', data=minVals)
    statsGroup.create_dataset("MAX_VALS", [Z], dtype='f4', data=maxVals)
    statsGroup.create_dataset("NAN_COUNTS", [Z], dtype='i4', data=nanCounts)

    histGroup = statsGroup.create_group("Histograms")
    histGroup.create_dataset("FIRST_CENTERS", [Z], dtype='f4', data=histogramFirstBinCenters)
    histGroup.create_dataset("BIN_WIDTHS", [Z], dtype='f4', data=histogramBinWidths)
    histGroup.create_dataset("BINS", [Z, N], dtype='i4', data=histogramBins)

    percentileGroup = statsGroup.create_group("Percentiles")
    percentileGroup.create_dataset("PERCENTILES", [len(percentiles)], dtype='f4', data=percentiles)
    percentileGroup.create_dataset("VALUES", [Z, len(percentiles)], dtype='f4', data=percentileVals)
    
    # WRITE AVERAGE_DATA_Z
    
    averageDataGroup = hduGroup.create_group("AverageData")
    dataSetAverageZ = averageDataGroup.create_dataset("AVERAGE_DATA_Z", dims[1:], dtype='f4', data=averageData)
    
    # WRITE AVERAGE_DATA_Z STATS
    
    averageZStatsGroup = allStatsGroup.create_group("StatisticsAverageZ")
    
    averageDataNaNFixed = averageData[~np.isnan(averageData)]
    means_avg_z = np.mean(averageDataNaNFixed)
    minVals_avg_z = np.min(averageDataNaNFixed)
    maxVals_avg_z = np.max(averageDataNaNFixed)
    
    nanCounts_avg_z = np.count_nonzero(np.isnan(averageData))
    
    percentileVals_avg_z = np.percentile(averageDataNaNFixed, percentiles)
    
    (tmpBins, tmpEdges) = np.histogram(averageDataNaNFixed, N)
    histogramBins_avg_z = tmpBins
    histogramBinWidths_avg_z = tmpEdges[1] - tmpEdges[0]
    histogramFirstBinCenters_avg_z = (tmpEdges[0] + tmpEdges[1]) / 2.0
    
    # TODO: this all needs to be factored out
    
    averageZStatsGroup.create_dataset("MEANS", [1], dtype='f4', data=means_avg_z)
    averageZStatsGroup.create_dataset("MIN_VALS", [1], dtype='f4', data=minVals_avg_z)
    averageZStatsGroup.create_dataset("MAX_VALS", [1], dtype='f4', data=maxVals_avg_z)
    averageZStatsGroup.create_dataset("NAN_COUNTS", [1], dtype='i4', data=nanCounts_avg_z)

    averageZhistGroup = averageZStatsGroup.create_group("Histograms")
    averageZhistGroup.create_dataset("FIRST_CENTERS", [1], dtype='f4', data=histogramFirstBinCenters_avg_z)
    averageZhistGroup.create_dataset("BIN_WIDTHS", [1], dtype='f4', data=histogramBinWidths_avg_z)
    averageZhistGroup.create_dataset("BINS", [1, N], dtype='i4', data=histogramBins_avg_z)

    averageZpercentileGroup = averageZStatsGroup.create_group("Percentiles")
    averageZpercentileGroup.create_dataset("PERCENTILES", [len(percentiles)], dtype='f4', data=percentiles)
    averageZpercentileGroup.create_dataset("VALUES", [1, len(percentiles)], dtype='f4', data=percentileVals_avg_z)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input filename')
    parser.add_argument('--chunks', nargs=3, type=int, help='Chunks to use, order: Z Y X')
    parser.add_argument('--stokes', help='Stokes parameter to assign', default='')
    parser.add_argument('--test', help='Run test and exit', action='store_true')
    args = parser.parse_args()
    
    if args.test:
        data = np.random.random(size=(5, 10, 15))
        data[data < 0.1] = np.nan
        data[0,:,:] = np.nan
        
        get_statistics(data)
        sys.exit()
    
    if args.stokes:
        print('Using Stokes parameter {}.'.format(args.stokes))
    else:
        print('No Stokes parameter specified.')
    
    baseFileName, _ = os.path.splitext(args.filename)

    if args.chunks:
        outputFileName = baseFileName + "_chunked_{}_{}_{}.hdf5".format(args.chunks[0], args.chunks[1], args.chunks[2])
    else:
        outputFileName = baseFileName + ".hdf5"
    
    with fits.open(args.filename) as inputFits:
        data = inputFits[0].data
        dims = data.shape
        
        if len(dims) == 4:
            dims = dims[1:]
            data = data[0, :, :, :]
            
        if len(dims) == 2:
            dims = (1, *dims)
            data = data[None, :]
            
        if len(dims) == 3:
            print('3D FITS file found, converting to HDF5 using IDIA customised LOFAR-USG-ICD-004 data structure')
            print('File dims: {}'.format(dims))
            print('Chunk dims: {}'.format(args.chunks))
            
            header = inputFits[0].header

            with h5py.File(outputFileName, "w") as outputHDF5:
                write_core(header, data, outputHDF5, args)
                write_statistics(header, data, outputHDF5, args)

        else:
            print('Only 3D FITS files supported for now')
            sys.exit(1)
