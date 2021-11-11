"""Install FEMDLpy and dependencies."""

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'A python package for a finite element based discretization of dynamic Laplacians based on the FEMDL MATLAB package'

setup(name='FEMDLpy',
    description=DESCRIPTION,
    author="Michael Denes",
    version=VERSION,
    packages=find_packages(),
    install_requires=[], #update
    keywords=['python', 'femdl', 'dynamic laplacian']
)
