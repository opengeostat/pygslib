"""
pygslib: GSLIB in python

Copyright 2015, Adrian Martinez Vargas.
Licensed under MIT.
"""

import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

from numpy.distutils.core import Extension

fgslib = Extension(name = 'pygslib.__fgslib',
                 sources = ['fgslib.f90'])

# This is a plug-in for setuptools that will invoke py.test
# when you run python setup.py test
class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest  
        sys.exit(pytest.main(self.test_args))

""" using this convention 
 major.minor[.build[.revision]]
 with development status at third position as follow: 
	0 for alpha (status)
	1 for beta (status)
	2 for release candidate
	3 for (final) release
"""
version = '0.0.0.2'

if __name__ == '__main__':
     
    #make sure you use the setup from numpy
	from numpy.distutils.core import setup
	setup(name='pygslib',
		  version=version,
		  description='Python wrap of GSLIB modified code and general geostatistical package',
		  long_description=open("README.rst").read(),
		  classifiers=[ 
		    'Development Status :: 3 - Alpha',
            'Programming Language :: Python',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: GIS'],
		  keywords='geostatistics kriging variogram estimation simulation', 
		  author='Adrian Martinez Vargas',
		  author_email='adriangeologo@yahoo.es',
		  url='https://github.com/opengeostat/pygslib',
		  license='MIT',
		  packages=find_packages(exclude=['examples', 'tests']),
		  include_package_data=True,
		  zip_safe=False,
		  tests_require=['numpy', 'pandas', 'matplotlib', 'nose'],
		  cmdclass={'test': PyTest},   
		  install_requires=['numpy', 'pandas', 'matplotlib'],
		  ext_modules = [fgslib])
