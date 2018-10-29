#!/bin/bash

conda config --add channels https://conda.anaconda.org/opengeostat
conda config --add channels conda-forge

pip install --no-deps ~/pygslib/pygslib_conda/wheels/pygslib-0.0.0.4.0.0-cp36-cp36m-linux_x86_64.whl
