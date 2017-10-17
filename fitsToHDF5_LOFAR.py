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


if len(sys.argv) < 2:
    print("Usage: {} <input file> <Stokes parameter to assign (opt)>".format(sys.argv[0]))
    quit(1)

inputFileName = sys.argv[1]
stokesParm = ''
if len(sys.argv) < 3:
    print('No Stokes parameter specified')
else:
    stokesParm = sys.argv[2][0]
    print('Using Stokes parameter {}'.format(stokesParm))

try:
    with fits.open(inputFileName) as inputFits:
        dims = inputFits[0].data.shape
        if len(dims) == 3:
            print('3D FITS file found, converting to HDF5 using LOFAR-USG-ICD-004 data structure')
            outputFileName = inputFileName + ".hdf5"
            outputFileName = outputFileName.replace(".fits", "")
            outputFileName = outputFileName.replace(".FITS", "")
            header = inputFits[0].header  # type: fits.Header

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
                        sliceData = inputFits[0].data[subBand, :, :]
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
                        polarizationCoordinateGroup.attrs.create('MultiStokes', False)
                        if len(stokesParm):
                            polarizationCoordinateGroup.attrs.create('StokesCoordinates', np.string_(stokesParm))

                        sourceTableGroup = currentGroup.create_group("SourceTable")
                        sourceHeaderAttributes = ['TELE', 'OBSERVER', 'INSTR', 'DATE-OBS', 'TIMESYS']
                        sourceHeaderAttributesFloat = ['OBSRA', 'OBSDEC', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z']
                        assignKeysToAttributes(sourceHeaderAttributes, header, sourceTableGroup)
                        assignKeysToAttributes(sourceHeaderAttributesFloat, header, sourceTableGroup, 'f')

                        processingHistoryGroup = currentGroup.create_group("ProcessingHistory")
                        dataSet = currentGroup.create_dataset("Data", sliceShape, dtype='f4', data=sliceData)
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
