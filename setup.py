# Install required Python version with Anaconda and Numpy
# conda create -n devdep python=3.10 numpy=2.1.1
# Removing Anaconda environmnet
# conda remove -n devdep --all
# To buil cython code in place use:
# python setup.py build_ext --inplace
# To run the main inside the module do
# python -m mymodule
# To buil the cython package:
# pip install .
# To remove gmapache:
# pip uninstall gmapache
# Add the following if running as a module
from distutils.core import setup
# Using setuptools and cython
from setuptools import setup
from Cython.Build import cythonize
# Turn on annotations with the following
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True
# Specify modules made with cython
setup(
    ext_modules = cythonize(
        ["gmapache/integerization.pyx"]
    )
)
