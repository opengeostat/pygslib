# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""PyGSLIB A module that wraps GSLIB fortran code into python.

This module provides an interface to FGSLIB (a python module 
generated automatically with f2p and modified gslib fortran code)

Copyright 2016, Adrian Martinez Vargas
                                                                        
This software may be modified and distributed under the terms        
of the MIT license.  See the LICENSE.txt file for details.  
         
"""                                                                      


import pandas as pd
import __gslib__kt3d
import __gslib__postik
import numpy as np
import matplotlib.pyplot as plt
import warnings

# Set default parameters in the Fortran common module

__gslib__kt3d.set_unest (np.nan) 
__gslib__postik.set_unest (np.nan)


#-----------------------------------------------------------------------------------------------------------------
#
#    Set nan values
#
#-----------------------------------------------------------------------------------------------------------------
def set_nan_value(value=np.nan):
    """Set non estimated value in the modules KT3D and POSTIK
    
    This will change the value assigned to non estimated values

    Parameters
    ----------
        value  :  numeric (optional, default np.nan) 
            
    Notes
    -----
    This will change the behaviour of KT3D in the actual python section.


    to see the actual nan_value you may call 
    
    __gslib__kt3d.UNEST
    
    """

    __gslib__kt3d.set_unest (value)
    __gslib__postik.set_unest (value)
    

#-----------------------------------------------------------------------------------------------------------------
#
#    Kriging with Kt3D
#
#-----------------------------------------------------------------------------------------------------------------
    
def kt3d(parameters):
    """Estimate with the GSLIB program KT3D
    
    This is a wrap for the GSLIB Fortran code of the program KT3D Version 2.0, 
    originally in Fortran 77. Only minor changes were included, the most 
    relevant are:
    
        a) the direct file output was redirected to Numpy Arrays
        b) the input (for grid estimate) is now as numpy arrays. 
           The grid definition was retained because is an input
           of the GSLIB search super-block algorithm.
        c) the trimming limits were removed, you may filter out 
           undesired values before estimating     
    
    the parameter file here is a dictionary with the following keys.
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, 
            some of the input are optional and will be internally initialized as
            zero or array of zeros. 

        kt3d_parameters = {
            # Input Data 
            # ----------
            'x' : ,   # rank-1 array('f') with bounds (nd)
            'y' : ,   # (optional) rank-1 array('f') with bounds (nd)
            'z' : ,   # (optional) rank-1 array('f') with bounds (nd)
            'vr' : ,   # (optional) rank-1 array('f') with bounds (nd)
            've' : ,   # (optional) rank-1 array('f') with bounds (nd)
            'bhid': , #  (optional) rank-1 array('i') with bounds (nd)
            # Output (Target) 
            # ----------
            'nx' : ,   # int
            'ny' : ,   # int
            'nz' : ,   # int
            'xmn' : ,   # float
            'ymn' : ,   # float
            'zmn' : ,   # float
            'xsiz' : ,   # float
            'ysiz' : ,   # float
            'zsiz' : ,   # float
            'nxdis' : ,   # int
            'nydis' : ,   # int
            'nzdis' : ,   # int
            'outx' : ,   # rank-1 array('f') with bounds (nout)
            'outy' : ,   # (optional) rank-1 array('f') with bounds (nout)
            'outz' : ,   # (optional) rank-1 array('f') with bounds (nout)
            'outextve' : ,   # (optional) rank-1 array('f') with bounds (nout)
            # Search parameters 
            # ----------
            'radius'     : ,   # float
            'radius1'    : ,   # float
            'radius2'    : ,   # float
            'ndmax'      : ,   # int
            'ndmin'      : ,   # int
            'noct'       : ,   # (optional) int
            'nbhid'       : ,   # (optional) int
            'sang1'      : ,   # (optional) float
            'sang2'      : ,   # (optional) float
            'sang3'      : ,   # (optional) float
            # Kriging parameters and options 
            # ----------
            'idrif'      : ,   # (optional) rank-1 array('i') with bounds (9)
            'itrend'     : ,   # (optional) int
            'ktype'      : ,   # (optional) int
            'skmean'     : ,   # (optional) float
            'koption'    : ,   # (optional) int
            'iktype'     : ,   # (optional) int
            'cut'        : ,   # (optional) rank-1 array('f') with bounds (ncut)
            'idbg'       : ,   # (optional) int
            # Variogram parameters 
            # ----------
            'c0'         : ,   # float
            'it'         : ,   # rank-1 array('i') with bounds (nst)
            'cc'         : ,   # rank-1 array('f') with bounds (nst)
            'aa'         : ,   # rank-1 array('f') with bounds (nst)
            'aa1'        : ,   # rank-1 array('f') with bounds (nst)
            'aa2'        : ,   # rank-1 array('f') with bounds (nst)
            'ang1'       : ,   # (optional) rank-1 array('f') with bounds (nst)
            'ang2'       : ,   # (optional) rank-1 array('f') with bounds (nst)
            'ang3'       : }   # (optional) rank-1 array('f') with bounds (nst)
            
        Returns
        -------
        dict
            estimation results at target points/blocks
            if iktype=0 returns variable estimate as:        
               {'outest'     : ,  #rank-1 array('f') with bounds (nout)
                'outkvar'    : }  #rank-1 array('f') with bounds (nout)
            if ktype=1 returns median kriging estimate as: 
               {'outcdf'     : }  #rank-2 array('f') with bounds (nout,ncut)
        dict
            debug output for the last block estimated     
            if idbg = 0
            {}
            if debug > 0 
            {'cbb'       : ,  #float
            'neq'        : ,  #int
            'na'         : ,  #int
            'dbgxdat'    : ,  #rank-1 array('f') with bounds (ndmax)
            'dbgydat'    : ,  #rank-1 array('f') with bounds (ndmax)
            'dbgzdat'    : ,  #rank-1 array('f') with bounds (ndmax)
            'dbgvrdat'   : ,  #rank-1 array('f') with bounds (ndmax)
            'dbgwt'      : ,  #rank-1 array('f') with bounds (ndmax + 11)
            'dbgxtg'     : ,  #float
            'dbgytg'     : ,  #float
            'dbgztg'     : ,  #float
            'dbgkvector' : ,  #rank-1 array('f') with bounds (ndmax + 11)
            'dbgkmatrix' : }  #rank-2 array('f') with bounds (ndmax + 11,ndmax + 11)
        dict
            estimation summary
            {} TODO: not implemented yet

      
    Notes
    -----
    Not all the parameters are used for calculation, it depend on the kriging type and options

    Optional parameters can be removed, knowing that kt3D will create an internal array/variable
    initialized to zero value
    
    If using nbhid > 0 the hole id number (bhid) is required. Hole IDs may be integers from
    one to total number of drillholes. Use function pygslib.Drillhole.txt2intID(table_name) to get 
    a correct bhid number.  

    This python functions is a wrapper that calls functions in the module __gslib__kt3d.so

    Example
    -------
    
    >>>     

    """

    # add dummy cut if not in parameters

    if 'cut' not in parameters:
        parameters['cut'] =[0]

    if 'cut' in parameters and parameters['cut']==None:
        parameters['cut'] =[0]
    
    # check that bhid is provided if nbhid > 0
    if parameters['nbhid']> 0 :
        assert parameters['bhid'] is not None, 'Error: BHID required if nbhid > 0'

    # check not using octants and min
    if parameters['nbhid']> 0 and parameters['noct']> 0:
        warnings.warn('/n Warning: !!!!! Using octants and maximum number of samples per drillholes at the same time /nmay produce unexpected results !!!!!!')

    # prepare the output
    output = {}
    debug  = {}
    estimate={}

    (output['outest'],     
    output['outkvar'],   
    output['outcdf'],     
    output['cbb'],        
    output['neq'],        
    output['na'],         
    output['dbgxdat'],    
    output['dbgydat'],    
    output['dbgzdat'],    
    output['dbgvrdat'],   
    output['dbgwt'],      
    output['dbgxtg'],    
    output['dbgytg'],     
    output['dbgztg'],     
    output['dbgkvector'], 
    output['dbgkmatrix'], error, warnings) =__gslib__kt3d.pykt3d(**parameters)



    if len(error.strip())>0:
        raise NameError(error.strip())

    if len(warnings.strip())>0:
        warnings.warn(warnings.strip())


    estimate['x']=parameters['outx']
    if 'outy' in parameters: 
        estimate['y']=parameters['outy']
    if 'outz' in parameters: 
        estimate['z']=parameters['outz']
    # clean a bit the output 
    if 'iktype' in parameters: 
        if parameters['iktype']==1:
            estimate['outcdf'] = output['outcdf']
        else: 
            estimate['outest']=output['outest']
            estimate['outkvar']=output['outkvar']
    else:
            estimate['outest']=output['outest']
            estimate['outkvar']=output['outkvar']
    
    debug ={}
    if 'idbg' in parameters:     
        if parameters['idbg']>0:
            for key in ['cbb','neq', 'na','dbgxtg', 'dbgytg', 'dbgztg']:
                debug[key]=output[key]
            debug['dbgxdat'] = output['dbgxdat'][:output['na']] 
            debug['dbgydat'] = output['dbgydat'][:output['na']]
            debug['dbgzdat'] = output['dbgzdat'][:output['na']]  
            debug['dbgvrdat'] = output['dbgvrdat'][:output['na']]
            debug['dbgwt'] = output['dbgwt'][:output['na']]
            debug['dbgkvector'] = output['dbgkvector'][:output['neq']]
            debug['dbgkmatrix'] = output['dbgkmatrix'][:output['neq'],:output['neq']]
            debug['lagrange'] = output['dbgwt'][output['na']:output['neq']]

            

            #plot data and weight

            fig, ax= plt.subplots()

            plt.axis((debug['dbgxtg'] - parameters['radius'],
                          debug['dbgxtg'] + parameters['radius'],
                          debug['dbgytg'] - parameters['radius'],
                          debug['dbgytg'] + parameters['radius']))

            plt.plot(parameters['x'],parameters['y'], '.', color = '0.25')    
            plt.scatter(debug['dbgxdat'], debug['dbgydat'], s=100+debug['dbgwt'], c=debug['dbgwt'], alpha=0.6)
            plt.plot (debug['dbgxtg'],debug['dbgytg'], 'r*')
            
            plt.colorbar()

            debug['plotxy'] = ax.get_figure()

    #TODO: Implement the estimation summary

    summary = {}

    return estimate, debug, summary




#-----------------------------------------------------------------------------------------------------------------
#
#    PostIK
#
#-----------------------------------------------------------------------------------------------------------------
    
def postik(parameters):
    """Post Process IK Distributions
    
    This is a wrap for the GSLIB Fortran code of the program POSTIK Version 2.0, 
    originally in Fortran 77. Only minor changes were included, the most 
    relevant are:
    
        a) the direct file output was redirected to Numpy Arrays
        b) the input (for grid estimate) is now as numpy arrays. 
        c) the trimming limits were removed, you may filter out 
           undesired values before post-processing    
    
    the parameter file here is a dictionary with the following keys.
     
    
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, 
            some of the input are optional and will be internally initialized as
            zero or array of zeros. 

        postik_parameters = {
            # output option, output parameter
            'iout'   : ,   # input int. Possible (1,2,3 or 4)
            'voutpar' : ,   # input float
            # the thresholds
            'ccut1'  : ,   # input rank-1 array('f') with bounds (nccut)
            # volume support?, type, varred
            'ivol'   : ,   # input int
            'ivtyp'  : ,   # input int
            'varred' : ,   # input float
            # minimum and maximum Z value
            'zmin'   : ,   # input float
            'zmax'   : ,   # input float
            # lower,middle and upper tail: option, parameter
            'ltail'  : ,   # input int
            'ltpar'  : ,   # input float
            'middle'  : ,   # input int
            'mpar'  : ,   # input float
            'utail'  : ,   # input int
            'utpar'  : ,   # input float
            # maximum discretization
            'maxdis' : ,   # input int
            # 1D arrays with global distribution
            'vr'     : ,   # input rank-1 array('f') with bounds (nc)
            'wt'     : ,   # input rank-1 array('f') with bounds (nc)
            # 2D array with IK3D output (continuous)
            'p'      : }   # input rank-2 array('f') with bounds (na,nc)
            
        Returns
        -------
            out1, out2, out3 : rank-1 array('f') with bounds (na) 
                results
                iout == 1 out1 -> 'mean', out2 -> NaN, out3 -> NaN
                iout == 2 out1 -> 'prob > cutoff', out2 -> 'mean > cutoff', out3 -> 'mean < cutoff'
                iout == 3 out1 -> 'Z value corresponding to CDF = value', out2 -> NaN, out3 -> NaN
                iout == 4 out1 -> 'Conditional Variance', out2 -> NaN, out3 -> NaN
               
      
    Notes
    -----
    iout options are: 
       1 = E-type
       2 = probability and mean above threshold(par)
       3 = Z percentile corresponding to (par)
       4 = conditional variance
    

    This python functions is a wrapper that calls functions in the module __gslib__postik.so

    Example
    -------
    
    >>>     out1,out2,out3,error = postik(parameters)

    """
  


    # prepare the output

    out1,out2,out3,error = __gslib__postik.postik(**parameters)

  
    if error>0:
        raise NameError('Error = {}'.format(error))

    return out1,out2,out3


