#!/bin/bash

conda config --add channels https://conda.anaconda.org/opengeostat

pip install --no-deps ~/pygslib/pygslib_conda/wheels/pygslib-0.0.0.4.0.0-cp27-cp27mu-linux_x86_64.whl
