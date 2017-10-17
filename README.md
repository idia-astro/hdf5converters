# hdf5converters
Python scripts to convert FITS and CASA images to HDF5 images using the LOFAR structure

## `fitsToHDF5_LOFAR.py`
`fitsToHDF5_LOFAR.py` converts a single 3D FITS file to the LOFAR-USG-ICD-004 data structure. Currently, the correct metadata is not written to file, but the file structure is created correctly, and an optional  `StokesCoordinates` attribute is written to the `PolarizationCoordinate` sub-group to identify which (if any) Stokes coordinate the dataset refers to. Metadata from the input FITS file is used to populate attributes of the relevant groups, but attribute names may differ from the LOFAR-USG-ICD-004 standard.

## `fitsToHDF5_LOFAR_4D.py`
`fitsToHDF5_LOFAR_4D.py` converts a list of 3D FITS files to a modified LOFAR-USG-ICD-004 data structure. Currently, the correct metadata is not written to file, but a `MultiStokes` flag and a list of the included Stokes coordinates in the form of the `StokesCoordinates` attribute are written to the `PolarizationCoordinate` sub-group to identify which Stokes coordinates the dataset includes. Metadata from the input FITS files is used to populate attributes of the relevant groups, but attribute names may differ from the LOFAR-USG-ICD-004 standard.
