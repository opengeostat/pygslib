"""
pygslib: GSLIB in python

Copyright 2020, Adrian Martinez Vargas.
Licensed under MIT.

"""
import setuptools  #required to build wheels
import sys
import os

# get version from package
exec(open('pygslib/version.py').read())

version = __version__
description = 'Python module for mineral resource estimation and geostatistics'
name='pygslib'
long_description=open("README.md").read()
classifiers=[
            'Development Status :: 6 - Alpha',
            'Programming Language :: Python',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License and GPL',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: GIS']
keywords='geostatistics kriging variogram estimation simulation blockmodel drillhole wireframe'
author='Adrian Martinez Vargas'
author_email='adriangeologo@yahoo.es'
url='https://github.com/opengeostat/pygslib'


if __name__ == '__main__':

    #FORTRAN code extension
    #-------------------------------------------------------------------
    #make sure you use the setup from numpy
    from numpy.distutils.core import setup # this is numpy's setup
    from numpy.distutils.core import Extension
    #from numpy import get_include
    from Cython.Build import cythonize

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
                     sources = ['for_code/rotscale.f90'],
								extra_link_args= ['-static'] )

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

    gslib_variograms = Extension(name = 'pygslib.gslib.__variograms',
                     sources = ['for_code/variograms.f90'] ) # extra_link_args=['-fbacktrace', '-fcheck=all']

    gslib_bigaus = Extension(name = 'pygslib.gslib.__bigaus',
                     sources = ['for_code/bigaus.f90'] )

    gslib_bicalib = Extension(name = 'pygslib.gslib.__bicalib',
                     sources = ['for_code/bicalib.f90'] )


    gslib_trans = Extension(name = 'pygslib.gslib.__trans',
                     sources = ['for_code/trans.f90'] )

    gslib_draw = Extension(name = 'pygslib.gslib.__draw',
                     sources = ['for_code/draw.f90'] )

    gslib_dm2csv = Extension(name = 'pygslib.gslib.__dm2csv',
                     sources = ['for_code/dm2csv.f90'] )


    # Cython

    drillhole = Extension(name ='pygslib.drillhole', sources =['cython_code/drillhole.pyx'])
    blockmodel= Extension(name ='pygslib.blockmodel', sources =['cython_code/blockmodel.pyx'])
    vtktools= Extension(name ='pygslib.vtktools', sources =['cython_code/vtktools.pyx'])
    nonlinear= Extension(name ='pygslib.nonlinear',sources =['cython_code/nonlinear.pyx'])
    sandbox= Extension(name ='pygslib.sandbox',sources =['cython_code/sandbox.pyx'])

    #pure python
    progress= Extension(name ='pygslib.progress',sources =['pygslib/progress.py'])
    surpac= Extension(name ='pygslib.surpac',sources =['pygslib/surpac.py'])

	# All extensions Fortran + Cython

    extensions =[gslib_variograms,
				 gslib_bigaus,
				 gslib_bicalib,
				 gslib_trans,
				 gslib_draw,
				 gslib_dm2csv,
				 gslib_addcoord,
				 gslib_rotscale,
				 gslib_read_gslib,
				 gslib_declus,
				 gslib_dist_transf,
				 gslib_block_covariance,
				 gslib_plot,
				 gslib_kt3d,
				 gslib_postik,
				 gslib_general,
				 drillhole,
				 blockmodel,
				 vtktools,
				 nonlinear,
				 sandbox,
                 progress,
				 surpac]

    extensions = cythonize(extensions, gdb_debug=True)

    setup(name=name,
          version=version,
          description= description,
          long_description=long_description,
          classifiers=classifiers,
          keywords=keywords,
          author=author,
          author_email=author_email,
          url=url,
          license='GPL and MIT',
          zip_safe=False,
          setup_requires = ['cython',
                            'numpy',
                            'setuptools',
                           ],
          install_requires = ['ipython',
                              'matplotlib',
                              'jupyter',
                              'vtk>=8.0',
                              'bokeh',
                              'colour',
                              'numpy>=0.19',
                              'scipy',
                              'pandas'
                              ],
          tests_require=['pytest'],
          packages=['pygslib',
                    'pygslib.gslib',
                    'pygslib.plothtml',
                    'pygslib.charttable'],
          include_package_data=True,
          package_data={'pygslib': ['data/*.*']},
          ext_modules = extensions)


	# copy dll if so is windows
    if os.name=='nt':
        mayor = sys.version_info[0]
        minor = sys.version_info[1]
        os.system('copy build\\lib.win-amd64-{}.{}\\pygslib\\.libs\\*.dll build\\lib.win-amd64-{}.{}\\pygslib\\gslib'.format(mayor,minor,mayor,minor))
        os.system('del /s /q build\\lib.win-amd64-{}.{}\\pygslib\\.libs'.format(mayor,minor))

    # running tests
    os.system('pytest --doctest-modules -v tests/') #pytest style tests

    print (" OPENGEOSTAT SAYS CYTHON/FORTRAN CODE COMPILED")
