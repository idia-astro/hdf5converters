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


if len(sys.argv) < 4:
    print("Usage: {} <stokes string> <input file 1> ... <input file N> <output file>".format(sys.argv[0]))
    quit(1)

stokesParms = sys.argv[1]
inputFileNames = sys.argv[2:-1]
if len(inputFileNames) != len(stokesParms):
    print("Number of input files ({}) does not match stokes string '{}'".format(len(inputFileNames), stokesParms))

inputFits = []
prevDims = []
for fileName in inputFileNames:
    try:
        nextInputFits = fits.open(fileName)
        dims = nextInputFits[0].data.shape
        if len(prevDims) and dims != prevDims:
            print('FITS files have inconsistent shapes!')
            quit(1)
        prevDims = dims
        inputFits.append(nextInputFits)

    except FileNotFoundError as error:
        print("File {} not found!".format(fileName))
        print(error)
        quit(1)

if len(prevDims) == 3:
    print('{} 3D FITS files found, converting to HDF5 using modified LOFAR-USG-ICD-004 data structure'.format(
        len(stokesParms)))
    outputFileName = sys.argv[-1]
    header = inputFits[0][0].header  # type: fits.Header
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

            imageGroupExpression = 'Image{0:03d}'
            if dims[0] > 999:
                imageGroupExpression = 'Image{0:04d}'
            for subBand in range(dims[0]):
                sliceShape = dims[1:]
                currentGroup = outputHDF5.create_group(imageGroupExpression.format(subBand))
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
                polarizationCoordinateGroup.attrs.create('MultiStokes', True)
                polarizationCoordinateGroup.attrs.create('StokesCoordinates', np.string_(stokesParms))

                sourceTableGroup = currentGroup.create_group("SourceTable")
                sourceHeaderAttributes = ['TELE', 'OBSERVER', 'INSTR', 'DATE-OBS', 'TIMESYS']
                sourceHeaderAttributesFloat = ['OBSRA', 'OBSDEC', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z']
                assignKeysToAttributes(sourceHeaderAttributes, header, sourceTableGroup)
                assignKeysToAttributes(sourceHeaderAttributesFloat, header, sourceTableGroup, 'f')

                processingHistoryGroup = currentGroup.create_group("ProcessingHistory")
                for stokesIndex in range(len(stokesParms)):
                    stokesVal = stokesParms[stokesIndex]
                    sliceData = inputFits[stokesIndex][0].data[subBand, :, :]
                    dataSet = currentGroup.create_dataset("Data{}".format(stokesIndex), sliceShape, dtype='f4', data=sliceData)
                    assignKeysToAttributes(['BUNIT'], header, dataSet)
                    dataSet.attrs.create('Stokes', np.string_(stokesVal))

                print("Finished processing {}".format(currentGroup.name))

    except OSError as error:
        print("Unable to create file {}".format(outputFileName))
        print(error)
