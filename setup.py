# to buil cython code in place use:
# python setup.py build_ext --inplace
# to buil cython inside a package just run:
# pip install .
# add the following if running as a module
from distutils.core import setup
# using setuptools and cython
from setuptools import setup
from Cython.Build import cythonize
# turn on annotations with the following
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True
# other packages, e.g. numpy
import numpy
# specify modules made with cython
setup(
    ext_modules = cythonize(
        ["gmapache/bla.pyx",
         "gmapache/ble.pyx",
         "gmapache/mysubmodule/bli.pyx"]
    ),
    include_dirs = [numpy.get_include()]
)
# when using numpy inside cython use the following and import numpy
# setup(
#     ext_modules = cythonize(
#         ["my_module.pyx"]
#     ),
#     include_dirs = [numpy.get_include()]
# )
