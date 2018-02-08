#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import sys
from typing import List
import argparse
import os


def convert(val, dtype=None):
    if dtype:
        return dtype(val)
    return np.string_(val)


def get_attr_copier(header):
    def copy_attrs(obj, keys: List[str], dtype=None):
        for key in keys:
            if key in header:
                obj.attrs.create(key, convert(header[key], dtype))
    return copy_attrs


def get_or_create_group(parent, name):
    if name in parent:
        return parent[name]
    return parent.create_group(name)


def write_core(inputFits, outputHDF5, args):
    header = inputFits[0].header
    copy_attrs = get_attr_copier(header)
    
    syslogGroup = outputHDF5.create_group("SysLog")
    
    if 'HISTORY' in header:
        syslogGroup.attrs.create('HISTORY', [np.string_(val) for val in header['HISTORY']])
        
    if 'COMMENT' in header:
        syslogGroup.attrs.create('COMMENT', [np.string_(val) for val in header['COMMENT']])

    currentGroup = get_or_create_group(outputHDF5, "Image")
    
    coordinatesGroup = currentGroup.create_group("Coordinates")
    
    directionCoordinatesGroup = coordinatesGroup.create_group("DirectionCoordinates")
    copy_attrs(directionCoordinatesGroup, ['CTYPE1', 'CUINIT1', 'CTYPE2', 'CUINIT2', 'RADESYS'])
    copy_attrs(directionCoordinatesGroup, ['CRVAL1', 'CRPIX1', 'CDELT1', 'CROTA1', 'CRVAL2', 'CRPIX2', 'CDELT2', 'CROTA2', 'EQUINOX'], float)

    spectralCoordinatesGroup = coordinatesGroup.create_group("SpectralCoordinate")
    copy_attrs(spectralCoordinatesGroup, ['CTYPE3', 'CUINIT3'])
    copy_attrs(spectralCoordinatesGroup, ['CRVAL3', 'CRPIX3', 'CDELT3', 'CROTA3'], float)

    polarizationCoordinateGroup = coordinatesGroup.create_group("PolarizationCoordinate")
    polarizationCoordinateGroup.attrs.create('MultiStokes', False)
    if args.stokes:
        polarizationCoordinateGroup.attrs.create('StokesCoordinates', convert(args.stokes))

    sourceTableGroup = currentGroup.create_group("SourceTable")
    copy_attrs(sourceTableGroup, ['TELE', 'OBSERVER', 'INSTR', 'DATE-OBS', 'TIMESYS'])
    copy_attrs(sourceTableGroup, ['OBSRA', 'OBSDEC', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z'], float)

    # Currently unused?
    #processingHistoryGroup = currentGroup.create_group("ProcessingHistory")
    
    if args.chunks:
        dataSet = currentGroup.create_dataset("Data", dims, dtype='f4', data=data, chunks=tuple(args.chunks))
    else:
        dataSet = currentGroup.create_dataset("Data", dims, dtype='f4', data=data)

    copy_attrs(dataSet, ['BUNIT'])

    if args.stokes:
        dataSet.attrs.create('Stokes', convert(args.stokes))


def write_statistics(inputFits, outputHDF5, args):
    data = inputFits[0].data
    dims = data.shape
    
    # CALCULATE

    percentiles = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 90.0, 95.0, 99.0, 99.5, 99.9, 99.95, 99.99, 99.999])

    means = np.zeros(dims[0] + 1)
    minVals = np.zeros(dims[0] + 1)
    maxVals = np.zeros(dims[0] + 1)
    nanCounts = np.zeros(dims[0] + 1)

    averageData = np.zeros(dims[1:])
    averageCount = np.zeros(dims[1:])

    percentileVals = np.zeros([dims[0] + 1, len(percentiles)])
    N = int(np.sqrt(dims[1] * dims[2]))
    histogramBins = np.zeros([dims[0] + 1, N])
    histogramFirstBinCenters = np.zeros(dims[0] + 1)
    histogramBinWidths = np.zeros(dims[0] + 1)

    for i in range(dims[0]):
        tmpData = data[i, :, :]
        nanArray = np.isnan(tmpData)
        nanCounts[i] = np.count_nonzero(nanArray)
        tmpDataNanFixed = tmpData[~nanArray]
        if nanCounts[i] == tmpData.shape[0] * tmpData.shape[1]:
            means[i] = np.NaN
            minVals[i] = np.NaN
            maxVals[i] = np.NaN
            percentileVals[i, :] = np.NaN
            histogramBins[i, :] = np.NaN
            histogramBinWidths[i] = np.NaN
            histogramFirstBinCenters[i] = np.NaN
        else:
            means[i] = np.mean(tmpDataNanFixed)
            minVals[i] = np.min(tmpDataNanFixed)
            maxVals[i] = np.max(tmpDataNanFixed)
            percentileVals[i, :] = np.percentile(tmpDataNanFixed, percentiles)
            (tmpBins, tmpEdges) = np.histogram(tmpDataNanFixed, N)
            histogramBins[i, :] = tmpBins
            histogramBinWidths[i] = tmpEdges[1] - tmpEdges[0]
            histogramFirstBinCenters[i] = (tmpEdges[0] + tmpEdges[1]) / 2.0

            averageCount += (~nanArray).astype(int)
            averageData += np.nan_to_num(tmpData)

    averageData /= np.fmax(averageCount, 1)
    averageData[averageCount < 1] = np.NaN
    averageDataNaNFixed = averageData[~np.isnan(averageData)]
    means[dims[0]] = np.mean(averageDataNaNFixed)
    minVals[dims[0]] = np.min(averageDataNaNFixed)
    maxVals[dims[0]] = np.max(averageDataNaNFixed)
    nanCounts[dims[0]] = np.count_nonzero(np.isnan(averageData))
    percentileVals[dims[0], :] = np.percentile(averageDataNaNFixed, percentiles)
    (tmpBins, tmpEdges) = np.histogram(averageDataNaNFixed, N)
    histogramBins[dims[0], :] = tmpBins
    histogramBinWidths[dims[0]] = tmpEdges[1] - tmpEdges[0]
    histogramFirstBinCenters[dims[0]] = (tmpEdges[0] + tmpEdges[1]) / 2.0
    
    # WRITE
    
    header = inputFits[0].header
    
    copy_attrs = get_attr_copier(header)

    statsGroup = outputHDF5.create_group("Statistics")

    statsGroup.create_dataset("Means", [dims[0] + 1], dtype='f4', data=means)
    statsGroup.create_dataset("MinVals", [dims[0] + 1], dtype='f4', data=minVals)
    statsGroup.create_dataset("MaxVals", [dims[0] + 1], dtype='f4', data=maxVals)
    statsGroup.create_dataset("NaNCounts", [dims[0] + 1], dtype='i4', data=nanCounts)

    histGroup = statsGroup.create_group("Histograms")
    histGroup.create_dataset("FirstCenters", [dims[0] + 1], dtype='f4', data=histogramFirstBinCenters)
    histGroup.create_dataset("BinWidths", [dims[0] + 1], dtype='f4', data=histogramBinWidths)
    histGroup.create_dataset("Bins", [dims[0] + 1, N], dtype='i4', data=histogramBins)

    percentileGroup = statsGroup.create_group("Percentiles")
    percentileGroup.create_dataset("Percentiles", [len(percentiles)], dtype='f4', data=percentiles)
    percentileGroup.create_dataset("Values", [dims[0] + 1, len(percentiles)], dtype='f4', data=percentileVals)
    
    currentGroup = get_or_create_group(outputHDF5, "Image")
        
    dataSetAverage = currentGroup.create_dataset("AverageData", dims[1:], dtype='f4', data=averageData)

    copy_attrs(dataSetAverage, ['BUNIT'])

    if args.stokes:
        dataSetAverage.attrs.create('Stokes', convert(args.stokes))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input filename')
    parser.add_argument('--chunks', nargs=3, type=int, help='Chunks to use, order: Z Y X')
    parser.add_argument('--stokes', help='Stokes parameter to assign', default='')
    args = parser.parse_args()
    
    if args.stokes:
        print('Using Stokes parameter {}'.format(args.stokes))
    else:
        print('No Stokes parameter specified')
    
    baseFileName, _ = os.path.splitext(args.filename)
    
    print(args.chunks)
    
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
            
        if len(dims) == 3:
            print('3D FITS file found, converting to HDF5 using IDIA customised LOFAR-USG-ICD-004 data structure')
            print('File dims: {}'.format(dims))
            print('Chunk dims: {}'.format(args.chunks))

            with h5py.File(outputFileName, "w") as outputHDF5:
                write_core(inputFits, outputHDF5, args)
                write_statistics(inputFits, outputHDF5, args)

        else:
            print('Only 3D FITS files supported for now')
            sys.exit(1)
