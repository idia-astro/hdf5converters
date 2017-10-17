from astropy.io import fits
import h5py
import numpy as np
import sys

if len(sys.argv)<2:
    print("Usage: {} <input file> <Stokes parameter to assign (opt)>".format(sys.argv[0]))
    quit(1)

inputFileName = sys.argv[1]
stokesParm = ''
if len(sys.argv)<3:
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

            try:
                with h5py.File(outputFileName, "w") as outputHDF5:
                    outputHDF5.create_group("SysLog")
                    imageGroupExpression = 'Image{0:03d}'
                    if dims[0] > 999:
                        imageGroupExpression = 'Image{0:04d}'
                    for subBand in range(dims[0]):
                        sliceShape = dims[1:]
                        sliceData = inputFits[0].data[subBand, :, :]
                        currentGroup = outputHDF5.create_group(imageGroupExpression.format(subBand))
                        coordinatesGroup = currentGroup.create_group("Coordinates")
                        directionCoordinatesGroup = coordinatesGroup.create_group("DirectionCoordinates")
                        spectralCoordinatesGroup = coordinatesGroup.create_group("SpectralCoordinate")
                        polarizationCoordinateGroup = coordinatesGroup.create_group("PolarizationCoordinate")
                        polarizationCoordinateGroup.attrs.create('MultiStokes', False)
                        if len(stokesParm):
                            polarizationCoordinateGroup.attrs.create('StokesCoordinates', np.string_(stokesParm))
                        sourceTableGroup = currentGroup.create_group("SourceTable")
                        processingHistoryGroup = currentGroup.create_group("ProcessingHistory")
                        dataSet = currentGroup.create_dataset("Data", sliceShape, dtype='f4', data=sliceData)
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
