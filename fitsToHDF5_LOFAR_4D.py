from astropy.io import fits
import h5py
import numpy as np
import sys

if len(sys.argv)< 4:
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
    print('{} 3D FITS files found, converting to HDF5 using modified LOFAR-USG-ICD-004 data structure'.format(len(stokesParms)))
    outputFileName = sys.argv[-1]

    try:
        with h5py.File(outputFileName, "w") as outputHDF5:
            outputHDF5.create_group("SysLog")
            imageGroupExpression = 'Image{0:03d}'
            if dims[0] > 999:
                imageGroupExpression = 'Image{0:04d}'
            #for subBand in range(dims[0]):
            for subBand in range(10):
                sliceShape = dims[1:]

                currentGroup = outputHDF5.create_group(imageGroupExpression.format(subBand))
                coordinatesGroup = currentGroup.create_group("Coordinates")
                directionCoordinatesGroup = coordinatesGroup.create_group("DirectionCoordinates")
                spectralCoordinatesGroup = coordinatesGroup.create_group("SpectralCoordinate")
                polarizationCoordinateGroup = coordinatesGroup.create_group("PolarizationCoordinate")
                sourceTableGroup = currentGroup.create_group("SourceTable")
                polarizationCoordinateGroup.attrs.create('MultiStokes', True)
                polarizationCoordinateGroup.attrs.create('StokesCoordinates', np.string_(stokesParms))
                processingHistoryGroup = currentGroup.create_group("ProcessingHistory")
                for stokesIndex in range(len(stokesParms)):
                    stokesVal = stokesParms[stokesIndex]
                    sliceData = inputFits[stokesIndex][0].data[subBand, :, :]
                    dataSet = currentGroup.create_dataset("Data{}".format(stokesIndex), sliceShape, dtype='f4', data=sliceData)
                print("Finished processing {}".format(currentGroup.name))


    except OSError as error:
        print("Unable to create file {}".format(outputFileName))
        print(error)

