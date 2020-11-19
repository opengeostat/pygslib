# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""PyGSLIB A module that wraps GSLIB fortran code into python.

This module provides an interface to GSLIB (a python module
generated automatically with f2p and modified gslib Fortran code)

Copyright 2015, Adrian Martinez Vargas

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE.txt file for details.

"""

# import 3d party modules
import sys
import os.path
import pandas as pd
import platform
import warnings
import numpy as np
import matplotlib.pyplot as plt

# import pygslib modules (compatible with python 3x)
from pygslib import drillhole, blockmodel, vtktools, nonlinear, sandbox, gslib, plothtml, charttable, progress, surpac

from pygslib import version

__version__= version.__version__


