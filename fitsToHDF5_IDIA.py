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
if chunks[0]<=0 or chunks[1]<=0 or chunks[2]<=0:
	doChunks = False

try:
	with fits.open(inputFileName) as inputFits:
		data = inputFits[0].data
		dims = data.shape
		if len(dims) == 4 and dims[0] ==1:
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

			try:
				with h5py.File(outputFileName, "w") as outputHDF5:
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
						#dataSetRotated = currentGroup.create_dataset("DataSwizzled", rotatedDims, dtype='f4', data=dataRotated, chunks=(64, 64, 16))
					else:
						dataSet = currentGroup.create_dataset("Data", dims, dtype='f4', data=data)
						#dataSetRotated = currentGroup.create_dataset("DataSwizzled", rotatedDims, dtype='f4', data=dataRotated)

					assignKeysToAttributes(['BUNIT'], header, dataSet)
					if len(stokesParm):
						dataSet.attrs.create('Stokes', np.string_(stokesParm))

						print("Finished processing {}".format(currentGroup.name))

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
