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

# Building the package and install dependencies inside anaconda environment
# pip install .

# Builing cython code in place if dependencies are already satisfied
# python setup.py build_ext --inplace
# Note - about Numpy Api Warning see:
# https://stackoverflow.com/questions/52749662/using-deprecated-numpy-api

# Running the main inside the module
# python -m gmapache

# Uninstalling gmapache from anaconda environment
# pip uninstall gmapache

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
    Extension("gmapache.partial_maps",
              ["gmapache/partial_maps.pyx"]),
    Extension("gmapache.integerization",
              ["gmapache/integerization.pyx"])
    # Extension("package.module",
    #           ["package/module.pyx"],
    #           include_dirs = [numpy.get_include()])
]

setup(
    name = "gmapache",
    ext_modules = cythonize(extensions)
)

################################################################################
################################################################################
