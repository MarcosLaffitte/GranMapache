################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications with                 #
#                Cython and HEuristics                                         #
#                                                                              #
# - Description: suite for the analysis of bijections, morphisms, alignments   #
#                and other maps between graphs.                                #
#                                                                              #
# - setup.py for gmapache package including cython                             #
#                                                                              #
################################################################################

# Installing required python version with anaconda
# conda create -n devdep python=3.11.10

# In order to install this package you need C and C++ compilers. If you dont have these
# already in your system you can install them INSIDE the anaconda environment with
# conda install -c conda-forge cxx-compiler

# Building the package and install dependencies inside anaconda environment
# python -m pip install .

# note that the last command is runned with "python -m pip", instead of just "pip",
# this allows to use the pip version linked to the required python version
# instead of the default version of the system

# Builing cython code in place if dependencies are already satisfied
# LINUX:   python setup.py build_ext --inplace
# WINDOWS: python setup.py build_ext -i --compiler=mingw32 -DMS_WIN64
# Note - about Numpy Api Warning see:
# https://stackoverflow.com/questions/52749662/using-deprecated-numpy-api

# Running the main inside the module
# python -m gmapache

# Uninstalling gmapache from anaconda environment
# python -m pip uninstall gmapache

# Removing anaconda environmnet completely
# conda remove -n devdep --all

# Using setuptools, cython and numpy
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

# Turn on annotations with the following
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

extensions = [
    Extension("gmapache.integerization",
              sources = ["gmapache/integerization.pyx"]),
    Extension("gmapache.partial_maps",
              sources = ["gmapache/partial_maps.pyx"],
              language = "c++",
              extra_compile_args=["-std=c++20"]),
    Extension("gmapache.isomorphisms",
              sources = ["gmapache/isomorphisms.pyx"],
              language = "c++",
              extra_compile_args = ["-std=c++20"]),
    # Extension("package.module",
    #           sources = ["package/module.pyx"],
    #           include_dirs = [numpy.get_include()],
    #           language = "c++",
    #           extra_compile_args=["-std=c++20"])
]

setup(
    name = "gmapache",
    ext_modules = cythonize(extensions)
)

################################################################################
################################################################################
