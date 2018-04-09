from setuptools import setup, find_packages

setup(
    name='idiahdf5',
    version='1.0',
    scripts=['scripts/fits2hdf5', 'scripts/fits2hdf5_4D', 'scripts/hdf5swizzle.py', 'scripts/hdf5swizzle_serial.py', 'scripts/profile_swizzle.sh'],
    packages=find_packages(),
)
