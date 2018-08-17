#!/bin/bash

conda config --add channels https://conda.anaconda.org/opengeostat

pip install --no-deps /home/adrian/pygslib/dist/pygslib-0.0.0.3.9-cp35-cp35m-linux_x86_64.whl

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
