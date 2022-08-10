# simple example from https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html
# compyle with python setup_cython.py build_ext --inplace

from numpy.distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("add_cython.pyx")
)