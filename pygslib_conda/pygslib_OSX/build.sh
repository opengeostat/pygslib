#!/bin/bash

conda config --add channels https://conda.anaconda.org/opengeostat

# this is to fix error in conda build... see https://github.com/numpy/numpy/issues/7427
LDFLAGS="$LDFLAGS -undefined dynamic_lookup -bundle"

$PYTHON setup.py install

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
