#!/bin/bash

conda config --add channels https://conda.anaconda.org/opengeostat
conda config --add channels conda-forge 

pip install --no-deps ~/Documents/OG_Python/pygslib/pygslib_conda/wheels/pygslib-0.0.0.6.0.0-cp38-cp38m-macosx_10_7_x86_64.whl
