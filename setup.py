# Install required Python version with Anaconda and Numpy
# conda create -n devdep python=3.11.10
# Removing Anaconda environmnet
# conda remove -n devdep --all
# To buil cython code in place use:
# python setup.py build_ext --inplace
# To run the main inside the module do
# python -m gmapache
# To buil the cython package:
# pip install .
# To remove gmapache:
# pip uninstall gmapache
# Using setuptools and cython
from setuptools import setup
from Cython.Build import cythonize
# Get numpy for usage inside cython
import numpy
# Turn on annotations with the following
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True
# Specify modules made with cython
setup(
    ext_modules = cythonize(
        ["gmapache/integerization.pyx"]
    ),
    include_dirs = [numpy.get_include()]
)
