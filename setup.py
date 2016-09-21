"""
pygslib: GSLIB in python

Copyright 2015, Adrian Martinez Vargas.
Licensed under MIT.

Note: 

In windows 64 bits use MinGW 64.

install dependencies: 

Microsoft Visual C++ Compiler for Python 2.7 at https://www.microsoft.com/en-ca/download/details.aspx?id=44266

c:\>conda install MinGW
c:\>conda install libpython

c:\>python setup.py config --compiler=mingw32 build --compiler=mingw32 install

"""
import warnings
import sys
from setuptools.command.test import test as TestCommand
from numpy.distutils.core import Extension


# check some dependencies

try:
 import vtk
except ImportError, e:
 warnings.warn('\nWarning:\n pygslib uses vtk but vtk is not installed!')


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

# define properties for setup
""" 
 Note: using this convention for version 
 major.minor[.build[.revision]]
 with development status at third position as follow: 
    0 for alpha (status)
    1 for beta (status)
    2 for release candidate
    3 for (final) release
"""
version = '0.0.0.3.8.4'
description = 'Python wrap of GSLIB modified code and general geostatistical package'
name='pygslib'
long_description=open("README.rst").read()
classifiers=[ 
            'Development Status :: 3 - Alpha',
            'Programming Language :: Python',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: GIS']
keywords='geostatistics kriging variogram estimation simulation'
author='Adrian Martinez Vargas'
author_email='adriangeologo@yahoo.es'
url='https://github.com/opengeostat/pygslib'
 
        
if __name__ == '__main__':
     
    #FORTRAN code extension
    #-------------------------------------------------------------------
    #make sure you use the setup from numpy
    from numpy.distutils.core import setup # this is numpy's setup
    from setuptools import find_packages

    # these are the almost intact gslib code
    gslib_kt3d = Extension(name = 'pygslib.gslib.__gslib__kt3d',
                     sources = ['for_code/kt3d/kt3d.f90', 
                                'for_code/kt3d/gslib/setrot.f90',
                                'for_code/kt3d/gslib/getindx.f90',
                                'for_code/kt3d/gslib/picksupr.f90',
                                'for_code/kt3d/gslib/setsupr.f90',
                                'for_code/kt3d/gslib/sqdist.f90',
                                'for_code/kt3d/gslib/cova3.f90',
                                'for_code/kt3d/gslib/ktsol.f90',
                                'for_code/kt3d/gslib/sortem.f90',
                                'for_code/kt3d/gslib/srchsupr.f90'],
                                f2py_options=[ 'only:', 'pykt3d', 'set_unest',  ':']) 

    gslib_postik = Extension(name = 'pygslib.gslib.__gslib__postik',
                     sources = ['for_code/postik/postik.f90', 
                                'for_code/postik/gslib/beyond.f90',
                                'for_code/postik/gslib/locate.f90',
                                'for_code/postik/gslib/powint.f90',
                                'for_code/postik/gslib/sortem.f90'],
                                f2py_options=[ 'only:', 'postik', 'set_unest', 'get_unest',  ':']) 
    
    
    # this is the gslib code too modified 
    # define extensions here:
    #-----------------------------------------------------  
    gslib_rotscale = Extension(name = 'pygslib.gslib.__rotscale',
                     sources = ['for_code/rotscale.f90'] )
                     
    gslib_block_covariance = Extension(name = 'pygslib.gslib.__block_covariance',
                     sources = ['for_code/block_covariance.f90'] )
                     
    gslib_read_gslib = Extension(name = 'pygslib.gslib.__read_gslib',
                     sources = ['for_code/read_gslib.f90'] )
                     
    gslib_addcoord = Extension(name = 'pygslib.gslib.__addcoord',
                     sources = ['for_code/addcoord.f90'] )                 
                     
    gslib_general = Extension(name = 'pygslib.gslib.__general',
                     sources = ['for_code/general.f90'] )  

    gslib_plot = Extension(name = 'pygslib.gslib.__plot',
                     sources = ['for_code/plot.f90'] ) 

    gslib_declus = Extension(name = 'pygslib.gslib.__declus',
                     sources = ['for_code/declus.f90'] ) 
                     
    gslib_dist_transf = Extension(name = 'pygslib.gslib.__dist_transf',
                     sources = ['for_code/dist_transf.f90'],
                     f2py_options=[ 'only:', 'backtr', 'anatbl',
                                    'nscore', 'ns_ttable', ':'] ) 
    # to exclude some fortran code use this: f2py_options=['only:', 'myfoo1', 'myfoo2', ':']  
    """
    dist_transf = Extension(name = 'pygslib.__dist_transf',
                     sources = ['for_code/dist_transf.f90'],
                     f2py_options=[ '--debug-capi',
                                    'only:', 'backtr', 
                                    'nscore', 'ns_ttable', ':'] ) 
    """

    gslib_variograms = Extension(name = 'pygslib.gslib.__variograms',
                     sources = ['for_code/variograms.f90'] ) 
                     
    gslib_bigaus = Extension(name = 'pygslib.gslib.__bigaus',
                     sources = ['for_code/bigaus.f90'] ) 

    gslib_bicalib = Extension(name = 'pygslib.gslib.__bicalib',
                     sources = ['for_code/bicalib.f90'] ) 

                     
    gslib_trans = Extension(name = 'pygslib.gslib.__trans',
                     sources = ['for_code/trans.f90'] ) 
                     
    gslib_draw = Extension(name = 'pygslib.gslib.__draw',
                     sources = ['for_code/draw.f90'] ) 
    
    # define fortran code setup in here 
    setup(name=name,
          version=version,
          description= description,
          long_description=long_description,
          classifiers=classifiers,
          keywords=keywords, 
          author=author,
          author_email=author_email,
          url=url,
          license='MIT',
          packages=find_packages(exclude=['examples', 'tests']),
          include_package_data=True,
          zip_safe=False,
          tests_require=['numpy', 'pandas', 'matplotlib', 'nose', 'mock'],
          cmdclass={'test': PyTest},   
          install_requires=['numpy', 'pandas', 'matplotlib'],
          ext_modules = [gslib_variograms,
                         gslib_bigaus,
                         gslib_bicalib,
                         gslib_trans,
                         gslib_draw,
                         gslib_addcoord,
                         gslib_rotscale,
                         gslib_read_gslib, 
                         gslib_declus,
                         gslib_dist_transf,
                         gslib_block_covariance,
                         gslib_plot,
                         gslib_kt3d,
                         gslib_postik,
                         gslib_general])
    
    print " OPENGEOSTAT SAYS FORTRAN CODE COMPILED"
	
    #Cython code extension
    #-------------------------------------------------------------------
    from distutils.core import setup    # this is the standard setup
    from distutils.extension import Extension as CYExtension
    #from Cython.Build import cythonize
    import numpy
    
    drillhole = CYExtension( 'pygslib.drillhole', 
                            ['cython_code/drillhole.pyx'], 
                            include_dirs=[numpy.get_include()]) 
                                
    blockmodel = CYExtension( 'pygslib.blockmodel', 
                            ['cython_code/blockmodel.pyx'],
                            include_dirs=[numpy.get_include()]) 
      
    vtktools = CYExtension( 'pygslib.vtktools', 
                            ['cython_code/vtktools.pyx'],
                            include_dirs=[numpy.get_include()]) 

    nonlinear = CYExtension( 'pygslib.nonlinear', 
                            ['cython_code/nonlinear.pyx'],
                            include_dirs=[numpy.get_include()]) 
    #deprecated
    #neighborhood = CYExtension( 'pygslib.neighborhood', 
                            #['cython_code/neighborhood.pyx'],
                            #include_dirs=[numpy.get_include()]) 

    #interpolators = CYExtension( 'pygslib.interpolators', 
                            #['cython_code/interpolators.pyx'],
                            #include_dirs=[numpy.get_include()])

                            
    setup(name=name,
          version=version,
          description= description,
          long_description=long_description,
          classifiers=classifiers,
          keywords=keywords, 
          author=author,
          author_email=author_email,
          url=url,
          license='GPL',
          packages=find_packages(exclude=['examples', 'tests']),
          include_package_data=True,
          zip_safe=False,
          tests_require=['numpy', 'pandas>=0.17', 'nose', 'mock', 'scipy'],
          cmdclass={'test': PyTest},   
          install_requires=['numpy', 'pandas', 'scipy'],
          ext_modules = [drillhole,
                        blockmodel,
                        vtktools,
                        nonlinear])

    print " OPENGEOSTAT SAYS CYTHON CODE COMPILED"