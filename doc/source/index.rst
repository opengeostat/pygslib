.. PyGSLIB Documentation documentation master file, created by
   sphinx-quickstart on Mon Aug 08 19:51:33 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyGSLIB documentation!
=================================================

PyGSLIB is an open source python package designed to do Mineral Resource Estimations with scripts.
It was inspired by Datamine Studio macros. Its philosophy is reproducibility and auditability.
With PyGSLIB you can write a script to do the entire resource estimate process,
from reading drillhole tables to estimation and validation.
You can rerun the entire script with different parameters until you get the desired output,
share it with your colleagues or send it to a peer reviewer for auditing.
PyGSLIB may be also useful to automatize an estimation process such as grade control model updates,
or to do loops required by some nonconventional applications such as drillhole
spacing studies with conditional simulations.


PyGSLIB is subdivided into five modules:

- ``gslib``. This is for geostatistics and interpolation. It was built with GSLIB FORTRAN code enhanced and linked to Python.
- ``drillhole``. This is for basic drillhole operation, such as compositing and desurveying.
- ``blockmodel``. This is for block modelling, it has functions to fill wireframes with blocks, reblocking, among others.
- ``vtktools``. This is for 3D computational geometry based on VTK, for example, to select samples within wireframes.  It also handles VTK files.
- ``nonlinear``. **This module is Under construction!**  It is an experimental module for nonlinear geostatistics based on the Discrete Gaussian Model.

Installation (Win64, Linux64, OSX64)
------------------------------------
The easiest way to install and work with PyGSLIB is using the Python distribution `Anaconda <https://www.continuum.io/downloads>`_.
To install PyGSLIB in the root environment of your anaconda distribution simply type in a terminal:

``conda install -c opengeostat pygslib``

You may also install `Paraview <http://www.paraview.org/download/>`_, this is required for visualization.

This video shows a demonstration of installation in Windows, these steps are similar for Linux or Mac installation.

.. raw:: html

    <div style="margin-top:10px;">
      <iframe width="560" height="315" src="https://www.youtube.com/embed/cbWXi7BfZVg" frameborder="0" allowfullscreen></iframe>
    </div>

Demonstration
--------------
You may read the :doc:`Tutorial  </Tutorial>` page and see the following video for a demonstration on resource estimation using PyGSLIB.


.. raw:: html

    <div style="margin-top:10px;">
      <iframe width="560" height="315" src="https://www.youtube.com/embed/SEwKy6wJbLE" frameborder="0" allowfullscreen></iframe>
    </div>

Contents:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Tutorial
   API

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
