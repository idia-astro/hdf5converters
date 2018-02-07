#!/usr/bin/env python3
from astropy.io import fits
import h5py
import numpy as np
import sys
from typing import List


def assignKeysToAttributes(keys: List[str], header: fits.Header, obj, dtype='str'):
	for key in keys:
		if key in header.keys():
			if dtype == 'f':
				obj.attrs.create(key, float(header[key]))
			elif dtype == 'i':
				obj.attrs.create(key, int(header[key]))
			else:
				obj.attrs.create(key, np.string_(header[key]))


if len(sys.argv) < 5:
	print("Usage: {} <input file> <Z-chunk> <Y-chunk> <X-chunk> <Stokes parameter to assign (opt)>".format(sys.argv[0]))
	quit(1)

inputFileName = sys.argv[1]
stokesParm = ''
if len(sys.argv) < 6:
	print('No Stokes parameter specified')
else:
	stokesParm = sys.argv[5][0]
	print('Using Stokes parameter {}'.format(stokesParm))

chunks = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
doChunks = True
if chunks[0] <= 0 or chunks[1] <= 0 or chunks[2] <= 0:
	doChunks = False

try:
	with fits.open(inputFileName) as inputFits:
		data = inputFits[0].data
		dims = data.shape
		if len(dims) == 4:
			dims = dims[1:]
			data = inputFits[0].data[0, :, :, :]
		if len(dims) == 3:
			print('3D FITS file found, converting to HDF5 using IDIA customised LOFAR-USG-ICD-004 data structure')
			print('File dims: {}'.format(dims))
			print('Chunk dims: {}'.format(chunks))
			if doChunks:
				outputFileName = inputFileName + "_chunked_{}_{}_{}.hdf5".format(chunks[0], chunks[1], chunks[2])
			else:
				outputFileName = inputFileName + ".hdf5"
			outputFileName = outputFileName.replace(".fits", "")
			outputFileName = outputFileName.replace(".FITS", "")
			header = inputFits[0].header  # type: fits.Header

			# rotatedDims = (dims[2], dims[1], dims[0])
			# dataRotated = np.zeros(rotatedDims, dtype='f4')
			# for i in range (rotatedDims[2]):
			# 	tmpData = np.array(inputFits[0].data[i, :, :])
			# 	dataRotated[:, :, i] = np.swapaxes(tmpData, 0, 1)
			# 	if i%10==0:
			# 		print('Rotated slice {}'.format(i))

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

			try:
				with h5py.File(outputFileName, "w") as outputHDF5:
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

					syslogGroup = outputHDF5.create_group("SysLog")
					if header['HISTORY']:
						historyVals = []
						for val in header['HISTORY']:
							historyVals.append(np.string_(val))
						syslogGroup.attrs.create('HISTORY', historyVals)

					if header['COMMENT']:
						commentVals = []
						for val in header['COMMENT']:
							commentVals.append(np.string_(val))
						syslogGroup.attrs.create('COMMENT', commentVals)

					currentGroup = outputHDF5.create_group("Image")
					coordinatesGroup = currentGroup.create_group("Coordinates")
					directionCoordinatesGroup = coordinatesGroup.create_group("DirectionCoordinates")
					assignKeysToAttributes(['CTYPE1', 'CUINIT1'], header, directionCoordinatesGroup)
					assignKeysToAttributes(['CRVAL1', 'CRPIX1', 'CDELT1', 'CROTA1'], header, directionCoordinatesGroup, 'f')
					assignKeysToAttributes(['CTYPE2', 'CUINIT2'], header, directionCoordinatesGroup)
					assignKeysToAttributes(['CRVAL2', 'CRPIX2', 'CDELT2', 'CROTA2'], header, directionCoordinatesGroup, 'f')
					assignKeysToAttributes(['RADESYS'], header, directionCoordinatesGroup)
					assignKeysToAttributes(['EQUINOX'], header, directionCoordinatesGroup, 'f')

					spectralCoordinatesGroup = coordinatesGroup.create_group("SpectralCoordinate")
					assignKeysToAttributes(['CTYPE3', 'CUINIT3'], header, spectralCoordinatesGroup)
					assignKeysToAttributes(['CRVAL3', 'CRPIX3', 'CDELT3', 'CROTA3'], header, spectralCoordinatesGroup, 'f')

					polarizationCoordinateGroup = coordinatesGroup.create_group("PolarizationCoordinate")
					polarizationCoordinateGroup.attrs.create('MultiStokes', False)
					if len(stokesParm):
						polarizationCoordinateGroup.attrs.create('StokesCoordinates', np.string_(stokesParm))

					sourceTableGroup = currentGroup.create_group("SourceTable")
					sourceHeaderAttributes = ['TELE', 'OBSERVER', 'INSTR', 'DATE-OBS', 'TIMESYS']
					sourceHeaderAttributesFloat = ['OBSRA', 'OBSDEC', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z']
					assignKeysToAttributes(sourceHeaderAttributes, header, sourceTableGroup)
					assignKeysToAttributes(sourceHeaderAttributesFloat, header, sourceTableGroup, 'f')

					processingHistoryGroup = currentGroup.create_group("ProcessingHistory")
					if doChunks:
						dataSet = currentGroup.create_dataset("Data", dims, dtype='f4', data=data, chunks=(64, 64, 16))
					# dataSetRotated = currentGroup.create_dataset("DataSwizzled", rotatedDims, dtype='f4', data=dataRotated, chunks=(64, 64, 16))
					else:
						dataSet = currentGroup.create_dataset("Data", dims, dtype='f4', data=data)
					# dataSetRotated = currentGroup.create_dataset("DataSwizzled", rotatedDims, dtype='f4', data=dataRotated)

					dataSetAverage = currentGroup.create_dataset("AverageData", dims[1:], dtype='f4', data=averageData)
					assignKeysToAttributes(['BUNIT'], header, dataSet)
					assignKeysToAttributes(['BUNIT'], header, dataSetAverage)
					if len(stokesParm):
						dataSet.attrs.create('Stokes', np.string_(stokesParm))
						dataSetAverage.attrs.create('Stokes', np.string_(stokesParm))


			except OSError as error:
				print("Unable to create file {}".format(outputFileName))
				print(error)
		else:
			print('Only 3D FITS files supported for now')
			quit(1)

except FileNotFoundError as error:
	print("File {} not found!".format(inputFileName))
	print(error)
	quit(1)
