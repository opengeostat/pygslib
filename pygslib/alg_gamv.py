# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""Algorithms for gamv


Copyright 2015, Adrian Martinez Vargas
                                                                        
This software may be modified and distributed under the terms        
of the MIT license.  See the LICENSE.txt file for details.  
         
"""                                                                      


#import os.path
#import pandas as pd
#import __fgslib
#import platform
#import warnings
import numpy as np
import matplotlib.pyplot as plt
import pygslib


#-----------------------------------------------------------------------------------------------------------------
#
#    Variograms GAMV algorithms
#
#-----------------------------------------------------------------------------------------------------------------

def vdirectional_1V(X,Y,Z,BHID,V1, tmin, tmax,
                            nlag,xlag,xltol,azm,atol,
                            bandwh,dip,dtol,bandwd,
                            isill,sills,ivtype):
    """Algorith to calculate univariate variogram in a given direction
    
    This function is to simplify the use of gamv for univariate variograms

    Parameters
    ----------
        X, Y , Z  :  1D array[float], 1D array[float], 1D array[float]
          Coordinate of in a cartesian space 
        BHID  :  1D array[int]
           Drillhole ID or Zone ID
        V1  :  1D array[float]
           Variable values to calculate the variogram
        tmin, tmax:  float, float
           trimming limits, float... default inf
        nlag,xlag,xltol: int, float, float
           number of lags, lag separation and lag tolerance 
           The default of xltol is xlag/2
        azm,atol,bandwh: float, float, float
           azimuth direction, azimuth tolerance and bandwith 'horizontal' 
           The default of bandwh is inf
        dip,dtol,bandwd: float, float, float
           dip, dip tolerance and bandwith 'vertical'
           The default of bandwd is inf
        isill, sills:  boolean, float  
           standarize sills?, variance used to std sills... 
           The default of sills is the variance of V1
       ivtype: int
           variogram type code

    Returns
    -------   
    out : ndarray, ndarray, ndarray,  ndarray, ndarray, matplotlib.fig
      The output is a set of 1D numpy arrays:
       pdis   Distance of pairs falling into this lag 
       pgam   Semivariogram, covariance, correlogram,... value
       phm    Mean of the tail data
       ptm    Mean of the head data
       pnump  Number of pairs
       fig    a matplotlib figure 

    Notes
    -----
    The output variogram type code are: 
    ivtype 1 = traditional semivariogram
           2 = traditional cross semivariogram
           3 = covariance
           4 = correlogram
           5 = general relative semivariogram
           6 = pairwise relative semivariogram
           7 = semivariogram of logarithms
           8 = semimadogram
    
    """
    
    #validate parameters and asign default values if necessary
    
    assert len(X)==len(Y)==len(Z)==len(BHID)==len(V1), "invalid array lengths"
    
    if xltol == None: 
        xltol= xlag/2.
    
    if bandwh == None:
        bandwh=np.inf
    
    if bandwd == None:
        bandwd=np.inf

    if isill==True and sills==None:
        sills=np.var(V1)

    if tmin == None:
        tmin=-np.inf
    
    if tmax == None:
        tmax=np.inf

    
    parameters = { 
    'x'      :  X,         # X coordinates, array('f') with bounds (nd), nd is number of data points
    'y'      :  Y,         # Y coordinates, array('f') with bounds (nd)
    'z'      :  Z,         # Z coordinates, array('f') with bounds (nd)
    'bhid'   :  BHID,      # bhid for downhole variogram, array('i') with bounds (nd)    
    'vr'     :  V1,        # Variables, array('f') with bounds (nd,nv), nv is number of variables
    'tmin'   :  tmin,      # trimming limits, float
    'tmax'   :  tmax,      # trimming limits, float
    'nlag'   :  nlag,      # number of lags, int
    'xlag'   :  xlag,      # lag separation distance, float                
    'xltol'  :  xltol,     # lag tolerance, float
    'azm'    :  [azm],       # azimut, array('f') with bounds (ndir)
    'atol'   :  [atol],      # azimut tolerance, array('f') with bounds (ndir)
    'bandwh' :  [bandwh],    # bandwith h, array('f') with bounds (ndir)
    'dip'    :  [dip],       # dip, array('f') with bounds (ndir)
    'dtol'   :  [dtol],      # dip tolerance, array('f') with bounds (ndir)
    'bandwd' :  [bandwd],    # bandwith d, array('f') with bounds (ndir)
    'isill'  :  isill,     # standardize sills? (0=no, 1=yes), int
    'sills'  :  [sills],     # variance used to std the sills, array('f') with bounds (nv)
    'ivtail' :  [1],       # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
    'ivhead' :  [1],       # head var., array('i') with bounds (nvarg)
    'ivtype' :  [ivtype],  # variogram type, array('i') with bounds (nvarg)
    'maxclp' :  1}         # maximum number of variogram point cloud to use, input int
        
    
    #test parameter files
    assert (pygslib.check_gamv_par(parameters))
    
    #calculate the variogram
    pdis,pgam, phm,ptm,phv,ptv,pnump, cldi, cldj, cldg, cldh=pygslib.gamv(parameters)
    
    #create plot
    v=0
    d=0
    width=xlag/2.
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()    
    ax2.bar(pdis[v, d, 1:]-width/2., pnump[v, d, 1:], width, zorder=-1)  
    ax2.set_ylabel('npoints', color='b') 
    ax2.set_zorder(-1)
    
    ax1.plot (pdis[v, d, 1:], pgam[v, d, 1:], '-o', label=str(dip) + '-->' + str(azm), color='r')
    ax1.set_xlabel('distance (h)')
    ax1.set_ylabel('gamma', color='r')

        
    return pdis[v, d, 1:], pgam [v, d, 1:], phm [v, d, 1:], ptm [v, d, 1:], pnump[v, d, 1:], fig

