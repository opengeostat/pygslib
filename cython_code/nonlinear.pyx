"""
PyGSLIB nonlinear, Module with function for nonlinear geostatistics  

Copyright (C) 2015 Adrian Martinez Vargas 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
any later version.
   
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Code based on paper

    A Step by Step Guide to Bi-Gaussian Disjunctive Kriging, by 
    Julian M. Ortiz, Bora Oz, Clayton V. Deutsch
    Geostatistics Banff 2004
    Volume 14 of the series Quantitative Geology and Geostatistics pp 1097-1102

See also
--------
 - http://www.ccgalberta.com/ccgresources/report05/2003-107-dk.pdf
 - http://www.ccgalberta.com/ccgresources/report04/2002-106-hermite.pdf
 - http://www.ccgalberta.com/ccgresources/report06/2004-112-inference_under_mg.pdf
 - Mining Geostatistics: A. G. Journel, Andre G. Journel, C. J
 - Introduction to disjunctive Kriging and non linear geostatistcs http://cg.ensmp.fr/bibliotheque/public/RIVOIRARD_Cours_00312.pdf

Warning
-------
This module is experimental, not tested and we are getting some validation
errors in the anamorphosis with declustering data. 

"""

cimport cython
cimport numpy as np
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
from libc.math cimport sqrt
from libc.math cimport exp
import pygslib
import matplotlib.pyplot as plt 


# is nan test for cython 
#from libc.math cimport isnan


# almost C version of stnormal 
@cython.boundscheck(False)
cdef float stnormpdf(float x):
    cdef float denom = (2*3.1415926)**.5
    cdef float num = exp(-x**2/2)
    return num/denom


#plotting options
ana_options = {
            'dt_pt_color' : 'grey',
            'dt_pt_line' : '-',
            'dt_pt_label' : 'data',
            'ex_pt_color' : 'black',
            'ex_pt_line' : '--',
            'ex_pt_label' : 'exp point',
            'ana_pt_color' : 'orange',
            'ana_pt_line' : '--',
            'ana_pt_label' : 'ana point',
            'ex_ana_pt_color' : 'green',
            'ex_ana_pt_line' : '-',
            'ex_ana_pt_label' : 'ana point(fixed)', 
            'ana_blk_color' : 'red',
            'ana_blk_line' : '--',
            'ana_blk_label' : 'ana block',
            'ex_ana_blk_color' : 'indigo',
            'ex_ana_blk_line' : '-', 
            'ex_ana_blk_label' : 'ana block(fixed)', }    
    
    
# ----------------------------------------------------------------------
#   Transformation table
# ----------------------------------------------------------------------
cpdef ttable(z, w = None):
    """ttable(z, w)
    
    Creates a transformation table. 
    
    Parameters
    ---------
    z : 1D numpy arrays of floats 
        Input variable.
    w : (Optional, default None) 1D numpy arrays of floats or None
        Variable and declustering weight. If None w is set internally as array of ones. 
    Returns:
    transin,transout : 1D numpy arrays with pairs of raw and gaussian values
    
    Note: 
    This function uses gslib.__dist_transf.ns_ttable
    """
    cdef int error
    
    if w is None:
        w = np.ones([len(z)])

    transin,transout, error = pygslib.gslib.__dist_transf.ns_ttable(z,w)
    
    assert error < 1, 'There was an error = {} in the function gslib.__dist_transf.ns_ttable'.format(error)
    
    
    return transin,transout
        
# ----------------------------------------------------------------------
#   Normal score
# ----------------------------------------------------------------------
cpdef nscore(z, transin, transout, getrank=False):
    """nscore(z, transin,transout, getrank)
    
    Normal score transformation, as in GSLIB
    
        
    Parameters
    ---------
    z : 1D numpy array of floats 
        Variable to transform.
    transin,transout : 1D numpy arrays of floats
        transformation table as obtained in function ttable
    getrank: boolean (optional, default False)
        If false return the gaussian values, if true returns gaussian probability (rank)

    Returns
    -------
    y : 1D numpy array of floats
        normal score transformation of variable z
    
    Note: 
    This function uses gslib.__dist_transf.nscore
    """   
    return pygslib.gslib.__dist_transf.nscore(z,transin,transout,getrank)
    
# ----------------------------------------------------------------------
#   Back transform using TTable 
# ----------------------------------------------------------------------
cpdef backtr(y, transin, transout, ltail, utail, ltpar, utpar, zmin, zmax, getrank=False):
    """nscore(z, transin,transout, getrank)
    
    Normal score transformation, as in GSLIB
    
        
    Parameters
    ---------
    y : 1D numpy array of floats 
        Gaussian values.
    transin,transout : 1D numpy arrays of floats
        transformation table as obtained in function ttable
    ltail : integer
    utail : integer
    ltpar : float
    utpar : float
    zmin : float
    zmax : float
    getrank : Boolean default False
    
    
    
    Returns
    -------
    z : 1D numpy array of floats
        raw transformation of gaussian variable y
    
    Note: 
    This function uses gslib.__dist_transf.backtr
    """   
    cdef int error
    
    z, error= pygslib.gslib.__dist_transf.backtr(vnsc = y, transin = transin, transout = transout, 
                                            ltail = ltail,
                                            utail = utail, # 4 is hyperbolic
                                            ltpar = ltpar,
                                            utpar = utpar,
                                            zmin = zmin,
                                            zmax = zmax, 
                                            getrank = getrank)
    
    assert error < 1, 'There was an error = {} in the function gslib.__dist_transf.backtr'.format(error)  

    return z
    
# ----------------------------------------------------------------------
#   Report some stas from a declusterd dataset
# ----------------------------------------------------------------------

# TODO: fix bug in quantiles reported
cpdef stats(z, w, iwt = True, report = True):
    """stats(z, w)
    
    Reports some basic stats using declustering weights 
    
    Parameters
    ---------
    z,w : 1D arrays of floats 
        Variable and declustering weight.
    iwt: boolean default True
        If True declustering weights will be used to calculate statistics
    report : boolean default True
        If True a printout will be produced

    Returns:
    xmin,xmax, xcvr,xmen,xvar : floats
        minimum, maximum, coefficient of variation, media and variance
    
    Note: 
    This function uses gslib.__plot.probplt to obtaining stats
    """
    cdef int error
    
    parameters_probplt = {
            'iwt'  : iwt,     #int, 1 use declustering weight
            'va'   : z,       # array('d') with bounds (nd)
            'wt'   : w}       # array('d') with bounds (nd), weight variable (obtained with declust?)



    binval,cl,xpt025,xlqt,xmed,xuqt,xpt975,xmin,xmax, xcvr,xmen,xvar,error = pygslib.gslib.__plot.probplt(**parameters_probplt)
    
    assert error < 1, 'There was an error = {} in the function gslib.__dist_transf.ns_ttable'.format(error)
    
    if report:
        print  ('Stats Summary')
        print  ('-------------')
        print  ('Count          ', len(z))
        print  ('Minimum        ', xmin)
        print  ('Maximum        ', xmax)
        print  ('CV             ', xcvr)
        print  ('Mean           ', xmen)
        print  ('Variance       ', xvar)
        print  ('Quantiles 2.5-97.5' , [xpt025,xlqt,xmed,xuqt,xpt975])
        
    return xmin,xmax, xcvr,xmen,xvar 
    
    
    
# ----------------------------------------------------------------------
#   Functions for punctual gaussian anamorphosis 
# ----------------------------------------------------------------------

#the recurrent formula for normalized polynomials
cpdef recurrentH(np.ndarray [double, ndim=1] y, int K=30):
    """recurrentH(np.ndarray [double, ndim=1] y, int K=30)
    
    Calculates the hermite polynomials with the recurrent formula
    
    Parameters
    ----------
    y : 1D array of float64
        Gaussian values calculated for the right part of the bin.
    K  : int32, default 30
        Number of hermite polynomials 

    Returns
    -------
    H : 2D array of float64
        Hermite monomials H(i,y) with shape [K+1,len(y)]
      
    See Also
    --------
    pygslib.gslib.__dist_transf.anatbl
       
    Note
    ----  
    The `y` values may be calculated on the right side of the bin, 
    as shown in fig VI.13, page 478 of Mining Geostatistics: 
    A. G. Journel, Andre G. Journel, C. J. The function 
    pygslib.gslib.__dist_transf.anatbl was prepared to provide these values,
    considering declustering weight if necessary. 
    
    The results from pygslib.gslib.__dist_transf.ns_ttable are inappropriate 
    for this calculation because are the mid interval in the bin.  
    """
    assert(K>=1)
    
    cdef np.ndarray [double, ndim=2] H
    cdef int k
    
    H=np.ones((K+1,len(y))) 
    #H[0,:]=1                #first monomial already ones 
    H[1,:]=-y               #second monomial
    
    # recurrent formula
    for k in range(1,K):
        H[k+1,:]= -1/np.sqrt(k+1)*y*H[k,:]-np.sqrt(k/float(k+1))*H[k-1,:]
    
    return H   #this is a 2D array of H (ki,yi)


#fit PCI for f(x)=Z(x)
cpdef fit_PCI(np.ndarray [double, ndim=1] z,
              np.ndarray [double, ndim=1] y,
              np.ndarray [double, ndim=2] H,
              double meanz=np.nan):
    """fit_PCI(np.ndarray [double, ndim=1] z, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=2] H, float meanz=np.nan)  
     
    Fits the hermite coefficient (PCI) 
    
    Parameters
    ----------
    z  : 1D array of float64
        Raw values sorted
    y  : 1D array of float64
        Gaussian values calculated for the right part of the bin.
    meanz: float64, default np.nan
        mean of z, if NaN then the mean will be calculated as np.mean(z)

    Returns
    -------
    PCI : 1D array of floats
        Hermite coefficients or PCI 
    g   : 1D array of floats
        pdf value (g[i]) corresponding to each gaussian value (y[i])
      
    See Also
    --------
    var_PCI
       
    Note
    ---- 
    ``PCI[0]=mean(z)`` and  ``sum=(PCI[1...n]^2)``. To validate the fit 
    calculate the variance with the function ``var_PCI()`` and compare 
    it with the experimental variance of `z`. You may also validate the 
    fit by calculating ``error= z-PHI(y)``, where ``PHI(y)`` are the 
    `z'` values calculated with the hermite polynomial expansion.  
    
    """
    
    assert y.shape[0]==z.shape[0]==H.shape[1], 'Error: wrong shape on input array'
    
    cdef np.ndarray [double, ndim=1] PCI
    cdef np.ndarray [double, ndim=1] g
    
    cdef unsigned int i, p, j, n=H.shape[0], m=H.shape[1]
    
    # if no mean provided
    if np.isnan(meanz):
        meanz = np.mean(z)
    
    PCI=np.zeros([H.shape[0]])
    g=np.zeros([H.shape[1]])
    PCI[0]=np.mean(z)
    
    for p in range(1,n):
        for i in range(1,m):
            g[i]= stnormpdf(y[i])
            PCI[p]=PCI[p] + (z[i-1]-z[i])*1/sqrt(p)*H[p-1,i]*g[i]
    
    return PCI, g


#get variance from PCI
cpdef var_PCI(np.ndarray [double, ndim=1] PCI):
    """var_PCI(np.ndarray [double, ndim=1] PCI) 
     
    Calculates the variance from hermite coefficient (PCI) 
     
    Parameters
    ----------
    PCI : 1D array of float64
        hermite coefficient

    Returns
    -------
    var : float64
        variance calculated with hermite polynomials
      
    See Also
    --------
    fit_PCI
       
    Note
    ----  
    The output may be used for validation of the PCI coefficients, it 
    may be close to the experimental variance of z.
    
    """
    
    a=PCI[1:]**2
    return np.sum(a)

#expand anamorphosis
cpdef expand_anamor(np.ndarray [double, ndim=1] PCI, 
                    np.ndarray [double, ndim=2] H,
                    double r=1.):
    """expand_anamor(np.ndarray [double, ndim=1] PCI, np.ndarray [double, ndim=2] H, double r=1.)
    
    Expands the anamorphosis function, that is :math:`Z = \sum_p(PSI_p*r^p*Hp(Yv))`
    
    r is the support effect. If r = 1 Z with point support will returned. 
    
    
    Parameters
    ----------
    PCI : 1D array of floats
        hermite coefficient
    H : 2D array of floats
        Hermite monomials H(i,y) with shape [K+1,len(y)]. See recurrentH
    r : float, default 1
        the support effect

    Returns
    -------
    PCI : 1D array of floats
        Hermite coefficients or PCI 
      
    See Also
    --------
    recurrentH
      
  
    """
    
    cdef np.ndarray [double, ndim=1] Z
    cdef int p
        
    Z=np.zeros(H.shape[1])
    Z[:]=PCI[0]
    for p in range(1,len(PCI)):
        Z+=PCI[p]*H[p,:]*r**p
    
    return Z

# ----------------------------------------------------------------------
#   Helper functions to preprocess punctual/experimental gaussian anamorphosis 
# ----------------------------------------------------------------------

     
# Back transformation from anamorphosis
cpdef Y2Z(np.ndarray [double, ndim=1] y,
        np.ndarray [double, ndim=1] PCI,
        double zamin, 
        double yamin, 
        double zpmin, 
        double ypmin, 
        double zpmax, 
        double ypmax, 
        double zamax, 
        double yamax,
        double r=1):
    """Y2Z(np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] PCI, double zamin, double yamin, double zpmin, double ypmin, double zpmax, double ypmax, double zamax, double yamax, double r=1)
    
    Gaussian (Y) to raw (Z) transformation 
    
    This is a convenience functions. It calls H=recurrentH(K,Y) and
    then returns Z = expand_anamor(PCI,H,r). K is deduced from 
    PCI.shape[0]. It also linearly interpolates the values
    out of the control points. 
    
    
    Parameters
    ----------
    PCI : 1D array of float64
        hermite coefficient
    y : 1D array of float64
        Gaussian values
    r : float64, default 1
        the support effect
    ypmin,zpmin,ypmax,zpmax : float64
         z, y practical minimum and maximum
    yamin,zamin,yamax,zamax : float64
         z, y authorized minimum and maximum

    Returns
    -------
    Z : 1D array of floats
        raw values corresponding to Y 
      
    See Also
    --------
    recurrentH, expand_anamor
       
    
    """
    
    cdef int K
    cdef np.ndarray [double, ndim=2] H
    cdef np.ndarray [double, ndim=1] Z
    cdef np.ndarray [double, ndim=1] zapmin= np.array([zamin,zpmin])
    cdef np.ndarray [double, ndim=1] yapmin= np.array([yamin,ypmin])
    cdef np.ndarray [double, ndim=1] zapmax= np.array([zpmax,zamax])
    cdef np.ndarray [double, ndim=1] yapmax= np.array([ypmax,yamax])

    
    K=PCI.shape[0]-1
    H=recurrentH(y,K)
    Z=expand_anamor(PCI,H,r)
    
    # fix some values based on the control points
    for i in range(y.shape[0]): 
        if y[i]<=ypmin:  
            Z[i]=np.interp(y[i], xp=yapmin, fp=zapmin)
            continue 
            
        if y[i]>=ypmax:  
            Z[i]=np.interp(y[i], xp=yapmax, fp=zapmax)
            continue 
        
    #and the new Z values with the existing PCI
    return Z

# Transformation from anamorphosis
cpdef Z2Y_linear(np.ndarray [double, ndim=1] z,
                 np.ndarray [double, ndim=1] zm,
                 np.ndarray [double, ndim=1] ym,
                 double zamin, 
                 double yamin, 
                 double zpmin, 
                 double ypmin, 
                 double zpmax, 
                 double ypmax, 
                 double zamax, 
                 double yamax):
    """Z2Y_linear(np.ndarray [double, ndim=1] z, np.ndarray [double, ndim=1] zm, np.ndarray [double, ndim=1] ym, double zamin, double yamin, double zpmin, double ypmin, double zpmax, double ypmax, double zamax, double yamax) 
             
    Raw (Z) to Gaussian (Y) transformation 
    
    Given a set of pairs [zm,ym] representing an experimental 
    Gaussian anamorphosis, this functions linearly interpolate y values 
    corresponding to z within the [zamin, zamax] intervals
    
    Parameters
    ----------
    z : 1D array of float64
        raw (Z) values where we want to know Gaussian (Y) equivalent
    zm,ym : 1D array of float64
        tabulated [Z,Y] values
    ypmin, zpmin, ypmax, zpmax : float64
         z, y practical minimum and maximum
    yamin, zamin, yamax,zamax : float64
         z, y authorized minimum and maximum

    Returns
    -------
    Y : 1D array of float64
        gaussian values corresponding to Z 
      
    See Also
    --------
    Y2Z
       
    
    """    
    
    cdef np.ndarray [double, ndim=1] Y=np.zeros(z.shape[0])
    cdef np.ndarray [double, ndim=1] zapmin= np.array([zamin,zpmin])
    cdef np.ndarray [double, ndim=1] yapmin= np.array([yamin,ypmin])
    cdef np.ndarray [double, ndim=1] zapmax= np.array([zpmax,zamax])
    cdef np.ndarray [double, ndim=1] yapmax= np.array([ypmax,yamax])

    
    # fix some values based on the control points
    for i in range(z.shape[0]): 
        
        if z[i]<=zamin:  
            Y[i]=yamin
            continue 

        if z[i]>=zamax:  
            Y[i]=yamax
            continue 
        
        if z[i]<=zpmin:  
            Y[i]=np.interp(z[i], xp=zapmin, fp=yapmin)
            continue 
            
        if z[i]>=zpmax:  
            Y[i]=np.interp(z[i], xp=zapmax, fp=yapmax)
            continue 
        
        if z[i]<zpmax and z[i]>zpmin:  
            Y[i]=np.interp(z[i], xp=zm, fp=ym)
            continue
        
    return Y

# ----------------------------------------------------------------------
#   Interactive gaussian anamorphosis modeling 
# ----------------------------------------------------------------------
cpdef calauthorized(zana, zraw, gauss, zpmin=None, zpmax=None):
    """calauthorized(zana, zraw, gauss, zpmin=None, zpmax=None)
    
    Calculate authorized intervals for gaussian anamorphosis with hermite polynomials
    
    Parameters
    ---------
    zana, zraw, gauss : 1D numpy arrays of floats 
        z values from anamorphosis, z values experimental and gaussian values.
        These arrays may have same shape and only one dimension. 
    zpmin, zpmax : floats (optional, default None)
        practical minimum and maximum values
    Returns:
    i, j, ii, jj : integer
        index of zamin, zamax, zpmin, zpmax on arrays zana, zraw, gauss
    
    Note: 
    This is a helper function used by pygslib.nonlinear.anamor
    """
    cdef int i
    cdef int j
    
    cdef int ii
    cdef int jj

    if zpmin is None:
        zpmin = min(zraw)

    if zpmax is None:
        zpmax = max(zraw)
        
    #get index for zpmax
    for jj in range(zraw.shape[0]-1, 0, -1):
        if zraw[jj]<zpmax:
            break
            
    #get index for zpmin
    for ii in range(0, zraw.shape[0]-1):
        if zraw[ii]>zpmin:
            break        
     
    # get index for minimum authorized
    for i in range(zana.shape[0]/2, 1, -1): 
        
        if zana[i-1] < zraw[ii] or zana[i-1] > zana[i] or gauss[i-1]<gauss[ii]:
            break
    
    # get index for maximum authorized
    for j in range(zana.shape[0]/2, zana.shape[0]-1, +1): 
        
        if zana[j+1] > zraw[jj] or zana[j+1] < zana[j] or gauss[j+1]>gauss[jj]:
            break 
        
    return i, j, ii, jj 

    
cpdef calauthorized_blk(zana, gauss, zpmin, zpmax):
    """calauthorized_blk(zana, gauss, zpmin, zpmax)
    
    Calculate authorized intervals for gaussian anamorphosis with hermite polynomials
    with block support. 
    
    Parameters
    ---------
    zana, gauss : 1D numpy arrays of floats 
        z values in block support from anamorphosis and gaussian values.
        These arrays may have same shape and only one dimension. 
    zpmin, zpmax : floats (optional, default None)
        practical minimum and maximum values
    Returns:
    i, j, ii, jj : integer
        index of zamin, zamax, zpmin, zpmax on arrays zana, zraw, gauss
    
    Note: 
    This is a helper function used by pygslib.nonlinear.anamor_blk.
    The index ii, jj of zpmin, zpmax are dum and set as jj = zana.shape[0]-1 and ii = 0
    
    """
    cdef int i
    cdef int j
    
    cdef int ii
    cdef int jj
        
    #get index for zpmax
    jj = zana.shape[0]-1
    #get index for zpmin
    ii = 0        
     
    # get index for minimum authorized
    for i in range(zana.shape[0]/2, 1, -1): 
        
        if zana[i-1] < zpmin or zana[i-1] > zana[i] or gauss[i-1]<gauss[ii]:
            break
    
    # get index for maximum authorized
    for j in range(zana.shape[0]/2, zana.shape[0]-1, +1): 
        
        if zana[j+1] > zpmax or zana[j+1] < zana[j] or gauss[j+1]>gauss[jj]:
            break 
        
    return i, j, ii, jj 

cpdef findcontrolpoints(zana, zraw, gauss, zpmin, zpmax, zamin, zamax):
    """findcontrolpoints(zana, zraw, gauss, zpmin, zpmax, zamin, zamax)
    
    Find the index of predefined authorized intervals in the arrays
    zana, zraw and gauss    
    
    Parameters
    ---------
    zana, zraw, gauss : 1D numpy arrays of floats 
        z values from anamorphosis, z values experimental and gaussian values.
        These arrays may have same shape and only one dimension.  
    zpmin, zpmax, zamin, zamax : floats
        practical minimum, practical maximum, authorized minimum and authorized maximum
    Returns:
    i, j, ii, jj : integer
        index of zamin, zamax, zpmin, zpmax on arrays zana, zraw, gauss
    
    Note: 
    This is a helper function used by pygslib.nonlinear.anamor.
   
    """
    cdef int i
    cdef int j

    cdef int ii
    cdef int jj
    
    assert zamax > zamin, 'zamax > zamin'
    assert zpmax > zpmin, 'zpmax > zpmin'
    assert zamax <= zpmax, 'zamax <= zpmax'
    assert zamin >= zpmin, 'zamin >= zpmin'
        
    #get index for zpmax
    for jj in range(zraw.shape[0]-1, 0, -1):
        if zraw[jj]<=zpmax:
            break
            
    #get index for zpmin
    for ii in range(0, zraw.shape[0]-1):
        if zraw[ii]>=zpmin:
            break        
     
    # get index for zamin
    for i in range(zana.shape[0]/2, 1, -1): 
        
        if zana[i-1] <= zamin:
            break
    
    # get index for zamax
    for j in range(zana.shape[0]/2, zana.shape[0]-1, +1): 
        
        if zana[j+1] >= zamax :
            break  

    assert zana[j] <= zraw[jj], 'Error: calculated zamax > calculated zpmax'
    assert zana[i] >= zraw[ii], 'Error: calculated zamin < calculated zpmin'
    assert gauss[j] <= gauss[jj], 'Error: calculated yamax > calculated ypmax'
    assert gauss[i] >= gauss[ii], 'Error: calculated yamin < calculated ypmin'    
    
    
    return i, j, ii, jj 


    
# Interactive anamorphosis modeling, including some plots 
def anamor(z, w, ltail=1, utail=1, ltpar=1, utpar=1, K=30, 
             zmin=None, zmax=None, ymin=None, ymax=None,
             zamin=None, zamax=None, zpmin=None, zpmax=None,
             ndisc = 1000, **kwargs):
    """anamor(z, w, ltail=1, utail=1, ltpar=1, utpar=1, K=30, 
             zmin=None, zmax=None, ymin=None, ymax=None,
             zamin=None, zamax=None, zpmin=None, zpmax=None,
             ndisc = 1000, **kwargs)
    
    Funtion to do interactive fitting of gaussian anamorphosis in point support as it follow

        - create experimental pairs z,y from input data
        - create an arbitrary secuence of Y ndisc values in interval [ymin, ymax] 
        - backtransform Y to obtain pairs Z, Y. You can define arbitrary 
          maximum and minimum practical values beyond data limits and extrapolation functions
        - calculate gaussian anamorphosis with hermite polynomials and hermite coefficients
        - calculate variance from hermite coefficients
        - asign or calculate the authorized and practical intervals of the gaussian anamorphosis
       
        You may change the parameters ltail, utail, ltpar, utpar, K, zmin, zmax, ymin, 
        ymax, zamin, zamax, zpmin, zpmax until the variance calculated with hermite coefficients 
        aproximates the declustered experimental variance.
        
        Note that maximum and minimum values zmin, zmax can be arbitrary and you use 
        zpmin, zpmax to create a transformation table within your desired intervals. 
    
    Parameters
    ---------
    z,w : 1D numpy arrays of floats 
        Variable and declustering weight.
    ltail, utail, ltpar, utpar: floats/integer (optional, default 1)
        extrapolation functions, as used in pygslib.nonlinear.backtr 
        and gslib.__dist_transf.backtr
    K: integer (optional, default 30)
        number of hermite polynomials.
    zmin, zmax: floats (optional, default None)
        minimum and maximum value in the lookup table used to fit the anamorphosis
        if None data limits are used, is these values are beyond data limits 
        interpolation functions are used to deduce Z values correponding to Y values
    ymin, ymax: floats (optional, default None)
        minimum and maximum gaussian values. If nont defined y equivalent of data 
        limits is used. Good values are -5 and 5, these are gaussian with low probability
    zamin, zamax, zpmin, zpmax: floats (optional, default None)
        Authorized and practical z value intervals. Authorized intervals are where 
        the gaussian anamorphosis do not fluctuate. Beyond authorized interval linear 
        extrapolation to practical interval are used. 
    ndisc: floats (optional, default 1000) 
        number of discretization intervals between ymin and ymax.
    **kwargs: more parameters
        extra parameters with ploting style
    
    Returns:
    PCI, H: hermite polynomial coefficients and monomials
    raw, zana, gauss: array of float
        z experimental extended, calculate with experimental anamorphosis and gaussian values
    Z: array of float
        Z of the gaussian anamorphosis pair [Z,Y] corrected
    P: array of float
        Probability P{Z<c}
    raw_var , PCI_var: float
        variance experimental (declustered) and calculated from hermite polynomials 
    fig: matplotlib figure
        gaussian anamorphosis plot 
        
    Note: 
    The pair [Z,P] defines the CDF in point support    
    
    """
    # set colors and line type
    options = ana_options
    options.update(kwargs)
    
    # z min and max for anamorphosis calculation
    if zmin==None:
        zmin = np.min(z)
        
    if zmax==None:
        zmax = np.max(z)
    
    # z min and max on any output... 
    if zpmin==None:
        zpmin = np.min(z)
        
    if zpmax==None:
        zpmax = np.max(z)
    
    
    # a) get experimental transformation table
    transin, transout = ttable(z, w)
    
    if ymin==None:
        ymin = np.min(transout)
        
    if ymax==None:
        ymax = np.max(transout)
    
    
    # b) generate a sequence of gaussian values
    gauss = np.linspace(ymin,ymax, ndisc)
    
    # c) Get the back transform using normal score transformation
    raw = pygslib.nonlinear.backtr(y = gauss, 
                                   transin = transin, 
                                   transout = transout, 
                                   ltail = ltail,
                                   utail = utail, # 4 is hyperbolic
                                   ltpar = ltpar,
                                   utpar = utpar,
                                   zmin = zmin,
                                   zmax = zmax,
                                   getrank = False)
                                   
    # d) get Hermite expansion
    H=  pygslib.nonlinear.recurrentH(y=gauss, K=K)
    
    # e) Get PCI
    xmin,xmax,xcvr,xmen,xvar = pygslib.nonlinear.stats(z,w,report=False)
    PCI, g =  pygslib.nonlinear.fit_PCI( raw, gauss,  H, meanz=np.nan)
    PCI[0] = xmen
    raw_var = xvar
    PCI_var = var_PCI(PCI)
    
    # f) Plot transformation
    zana = expand_anamor(PCI, H, r=1)
    
    if zamax is None or zamin is None:
        i, j, ii, jj = calauthorized(zana=zana,
                                     zraw=raw, 
                                     gauss=gauss,
                                     zpmin=zpmin, 
                                     zpmax=zpmax)
    else: 
        i, j, ii, jj = findcontrolpoints(zana=zana,
                                     zraw  = raw, 
                                     gauss = gauss,
                                     zpmin = zpmin, 
                                     zpmax = zpmax,
                                     zamin = zamin,
                                     zamax = zamax)

    # d) Now we get the transformation table for gaussian anamorphosis corrected
    Z = pygslib.nonlinear.backtr(  y = gauss, 
                                   transin = zana[i:j+1], 
                                   transout = gauss[i:j+1], 
                                   ltail = ltail,
                                   utail = utail, # 4 is hyperbolic
                                   ltpar = ltpar,
                                   utpar = utpar,
                                   zmin = raw[ii],
                                   zmax = raw[jj],
                                   getrank = False)
    
    # e) And the probability  Prob[Z<cz] knowing that Prob[Z<cz] == Prob[Y<cz]
    P = norm.cdf(gauss)
    
    # plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylabel('Z')
    plt.xlabel('Y')
    ax.plot(transout, transin, options['dt_pt_line'], color = options['dt_pt_color'], linewidth=4.0, label = options['dt_pt_label'])
    ax.plot(gauss[1:-1],raw[1:-1], options['ex_pt_line'], color = options['ex_pt_color'], label = options['ex_pt_label'])
    ax.plot(gauss[1:-1],zana[1:-1], options['ana_pt_line'], color = options['ana_pt_color'], label = options['ana_pt_label'])
    ax.plot(gauss,Z, options['ex_ana_pt_line'], color = options['ex_ana_pt_color'], label = options['ex_ana_pt_label'])
    ax.plot(gauss[i],zana[i], 'or',mfc='none')
    ax.plot(gauss[j],zana[j], 'or',mfc='none')
    ax.plot(gauss[ii],raw[ii], 'ob',mfc='none')
    ax.plot(gauss[jj],raw[jj], 'ob',mfc='none')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    
    print ('Raw Variance', raw_var)
    print ('Variance from PCI', PCI_var)
    
    print ('zamin', zana[i])
    print ('zamax', zana[j])
    print ('yamin', gauss[i])
    print ('yamax', gauss[j])
    
    print ('zpmin', raw[ii])
    print ('zpmax', raw[jj])
    print ('ypmin', gauss[ii])
    print ('ypmax', gauss[jj])
    
    
    return PCI, H, raw, zana, gauss,Z , P, raw_var , PCI_var, fig

# interactive block anamorphosis transformation
def anamor_blk( PCI, H, r, gauss, Z,
                  ltail=1, utail=1, ltpar=1, utpar=1,
                  raw=None, zana=None, **kwargs):  
    """anamor_blk( PCI, H, r, gauss, Z,
                  ltail=1, utail=1, ltpar=1, utpar=1,
                  raw=None, zana=None, **kwargs)
    
    Funtion to do obtain gaussian anamorphosis in block support from punctual anamorphosis
     
    
    Parameters
    ---------
    PCI, H: 
        hermite polynomial coefficients and monomials
    r: float
        support correlation coefficient
    gauss, Z: array of floats
        pair of gaussian and Z obtained with anamorphosis
    ltail, utail, ltpar, utpar: floats/integer (optional, default 1)
        extrapolation functions, as used in pygslib.nonlinear.backtr 
        and gslib.__dist_transf.backtr
    raw, zana: array of floats (default None)
        extended experimental and non-transformed gaussian z calculated with hermite polynomials
    **kwargs: more parameters
        extra parameters with ploting style
    
    Returns:

    ZV: array of float
        corrected Z values in block support
    PV: array of float
        Probability P{Z<c}
    fig: matplotlib figure
        gaussian anamorphosis plot     
    
    Note: 
    The pair [ZV,PV] defines the CDF in block support
        
    """
    
    # set colors and line type
    options = ana_options
    options.update(kwargs)
    
    # get practical limits on Z
    zpmin = np.min(Z)
    zpmax = np.max(Z)
    
    # Get Z experimental
    z_v= expand_anamor(PCI, H, r)
    
    # Get authorized limits on z experimental
    i, j, ii, jj = calauthorized_blk( zana=z_v, 
                                     gauss=gauss,
                                     zpmin=zpmin, 
                                     zpmax=zpmax)

    

    print ('zamin blk', Z[i])
    print ('zamax blk', Z[j])
    print ('yamin blk', gauss[i])
    print ('yamax blk', gauss[j])
    
    print ('zpmin blk', zpmin)
    print ('zpmax blk', zpmax)
    print ('ypmin blk', gauss[ii])
    print ('ypmax blk', gauss[jj])


    # Now we get the transformation table corrected
    ZV = pygslib.nonlinear.backtr( y = gauss, 
                                   transin = z_v[i:j+1], 
                                   transout = gauss[i:j+1], 
                                   ltail = ltail,
                                   utail = utail, # 4 is hyperbolic
                                   ltpar = ltpar,
                                   utpar = utpar,
                                   zmin = zpmin,
                                   zmax = zpmax,
                                   getrank = False)

    # And the probability  Prob[ZV<cz] knowing that Prob[ZV<cz] == Prob[YV<cz]
    PV = norm.cdf(gauss)                                                            
    
    # plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylabel('Z')
    plt.xlabel('Y')
    
    
    
    if raw is not None: 
        ax.plot(gauss[1:-1],raw[1:-1], options['ex_pt_line'], color = options['ex_pt_color'], label =  options['ex_pt_label'])
    if zana is not None:
        ax.plot(gauss[1:-1],zana[1:-1], options['ana_pt_line'], color = options['ana_pt_color'], label =  options['ana_pt_label'])
    
    ax.plot(gauss[1:-1],Z[1:-1], options['ex_ana_pt_line'], color = options['ex_ana_pt_color'], label =  options['ex_ana_pt_label'])
    
    ax.plot(gauss[1:-1],z_v[1:-1], options['ana_blk_line'], color = options['ana_blk_color'], label =  options['ana_blk_label'])
    
    ax.plot(gauss,ZV, options['ex_ana_blk_line'], color = options['ex_ana_blk_color'], label =  options['ex_ana_blk_label'])
    
    ax.plot(gauss[i],z_v[i], 'or',mfc='none')
    ax.plot(gauss[j],z_v[j], 'or',mfc='none')
    ax.plot(gauss[ii],ZV[ii], 'ob',mfc='none')
    ax.plot(gauss[jj],ZV[jj], 'ob',mfc='none')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    return ZV,PV, fig
    
    
# Direct anamorphosis modeling from raw data
def anamor_raw(z, w, K=30, **kwargs):
    """anamor_raw(z, w, K=30, **kwargs)
    
    Non-interactive gaussian anamorphosis calculated directly from data 
        
    Parameters
    ---------
    z,w : 1D numpy arrays of floats 
        Variable and declustering weight.
    K: integer (optional, default 30)
        number of hermite polynomials.
    **kwargs: more parameters
        extra parameters with ploting style
    
    Returns:
    PCI, H: hermite polynomial coefficients and monomials
    raw, zana, gauss: array of float
        z experimental extended, calculate with experimental anamorphosis and gaussian values
    raw_var , PCI_var: float
        variance experimental (declustered) and calculated from hermite polynomials 
    fig: matplotlib figure
        gaussian anamorphosis plot 
    
    Note: 
    
    
    """ 
    # set colors and line type
    options = ana_options
    options.update(kwargs)

    
    # a) get experimental transformation table
    raw, gauss = ttable(z, w)   
                                   
    # b) get Hermite expansion
    H=  pygslib.nonlinear.recurrentH(y=gauss, K=K)
    
    # e) Get PCI
    xmin,xmax,xcvr,xmen,xvar = pygslib.nonlinear.stats(z,w,report=False)
    PCI, g =  pygslib.nonlinear.fit_PCI( raw, gauss,  H, meanz=np.nan)
    PCI[0] = xmen
    raw_var = xvar
    PCI_var = var_PCI(PCI)
    
    # f) Plot transformation
    zana = pygslib.nonlinear.Y2Z(gauss, PCI, 
                           zamin = raw.min(), 
                           yamin = gauss.min(),  
                           zpmin = raw.min(),  
                           ypmin = gauss.min(), 
                           zpmax = raw.max(), 
                           ypmax = gauss.max(), 
                           zamax = raw.max(),
                           yamax = gauss.max(),
                           r=1.)    

                           
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylabel('Z')
    plt.xlabel('Y')
    ax.plot(gauss,raw, options['dt_pt_line'], color = options['dt_pt_color'], linewidth=4.0, label = options['dt_pt_label'])
    ax.plot(gauss,zana, options['ana_pt_line'], color = options['ana_pt_color'], label = options['ana_pt_label'])
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    print ('Raw Variance', raw_var)
    print ('Variance from PCI', PCI_var)
    
    return PCI, H, raw, zana, gauss, raw_var , PCI_var, fig
    

# ----------------------------------------------------------------------
#   Extra Functions for support and information effect  
# ----------------------------------------------------------------------
cpdef f_var_Zv(double r,
               np.ndarray [double, ndim=1] PCI,
               double Var_Zv=0):
    """f_var_Zv(double r, np.ndarray [double, ndim=1] PCI, double Var_Zv=0)
    
    This is an internal function used to deduce the coefficients r
    (or s), representing the support effect. It defines the relations:  
    
    
        :math:`Var(Zv) = \sum PCI^2 * r^(n*2)`
        
        or 
    
        :math:`Var(Zv*) = \sum PCI^2 * s^(n*2)`
    
    r is necessary to account for information effect
    s is necessary to account for smoothing in the information effect.        
        
    see "The information effect and estimating recoverable reserves"
    J. Deraisme (Geovariances), C. Roth (Geovariances) for more information
    
    Parameters
    ----------
    r : float64
        r or s coefficient representing support effect of Zv (or Zv*)
    PCI : 1D array of floats
        hermite coefficients
    Var_Zv : float64
        Block Variance var(Zv) or var(Zv*) 
    
    Note
    ----
    var(Zv) can be calculated as C(0)-gamma(v,v) or C(v,v) see function 
    block_covariance 
    
    var(Zv*) = var(Zv) - Kriging variance - 2*LaGrange multiplier
     
    In case of information effect this can be calculated with a dummy 
    dataset in a single block representing future information, for 
    example blast holes. 
    
    
    """
    
    cdef float a 
    cdef int i
    
    a=0.
    
    for i in range(1,len(PCI)):
       a+=PCI[i]**2. * r**(2.*i)
    
    
    return a-Var_Zv

# this is to pass a function to brentq
# auxiliar function covar (Zv,Zv*) = sum PCI^2 * r^n * s^n * ro^n 
# see "The information effect and estimating recoverable reserves"
# J. Deraisme (Geovariances), C. Roth (Geovariances)
cpdef f_covar_ZvZv(double ro,
                   double s,
                   double r,
                   np.ndarray [double, ndim=1] PCI,
                   double Covar_ZvZv=0):
    """f_covar_ZvZv(double ro, double s, double r, np.ndarray [double, ndim=1] PCI, double Covar_ZvZv=0)
    
    This is an internal function used to deduce the coefficients 
    ro = covar(Yv, Yv*). This function represents the expression:  
    
    
        :math:`Covar (Zv,Zv^*) = \sum PCI^2 * r^n * s^n * ro^n`
        
    ro is necessary to account for the conditional bias in the 
    information effect.     
        
    see "The information effect and estimating recoverable reserves"
    J. Deraisme (Geovariances), C. Roth (Geovariances) for more information
    
    Parameters
    ----------
    r, ro, s : float
        support effect and information effect coefficients.
    PCI : 1D array of floats
        hermite coefficients
    Covar_ZvZv : float
        Block covariance (correlation) between true Zv and estimate Zv* 
    
    Note
    ----
    :math:`Covar (Zv,Zv^*) = var(Zv) - Kriging variance - LaGrange multiplier`
    
    see expression 7.47, page 97 on Basic Linear Geostatistics by 
    Margaret Armstrong.
    
    In case of information effect this can be calculated with a dummy 
    dataset in a single block representing future information, for 
    example blast holes. 
    
    Note that the slop of regression is  
    
    :math:`p = Covar (Zv,Zv*) / (Covar (Zv,Zv*)  - LaGrange multiplier)`
    
    """
    
    cdef float a 
    cdef int i
    
    a=0.
    
    for i in range(1,len(PCI)):
       a+=PCI[i]**2. * r**i * s**i * ro**i
    
    
    return a-Covar_ZvZv


#calculate support effect coefficient r
cpdef get_r(double Var_Zv,
            np.ndarray [double, ndim=1] PCI):
    """get_r(double Var_Zv, np.ndarray [double, ndim=1] PCI)
    
    This function deduces the value of the support effect coefficient r
    or the information effect coefficient, smoothing component, s 
    defined by the equations: 
    
        :math:`Var(Zv) = \sum PCI^2 * r^(n*2)`
        
        and
    
        :math:`Var(Zv*) = \sum PCI^2 * s^(n*2)`
    
    
    The value of r is deduced by finding the root of the equation 
    f_var_Zv, using the classic Brent method (see scipy.optimize.brentq) 
    
    
    Parameters
    ----------
    PCI : 1D array of float64
        hermite coefficient
    Var_Zv : float64
        Block variance

    Returns
    -------
    r :  float64
        support effect coefficient r or information effect coefficient s
      
    See Also
    --------
    f_var_Zv, fit_PCI, anamor, scipy.optimize.brentq
    
    Note
    ----
    var(Zv) can be calculated as C(0)-gamma(v,v) or C(v,v) see function 
    block_covariance 
    
    var(Zv*) = var(Zv) - Kriging variance - 2*LaGrange multiplier
     
    In case of information effect this can be calculated with a dummy 
    dataset in a single block representing future information, for 
    example blast holes. 
    
    """
    
    return brentq(f=f_var_Zv, a=0, b=1, args=(PCI,Var_Zv))

#calculate information effect coefficient ro
cpdef get_ro(double Covar_ZvZv,
           np.ndarray [double, ndim=1] PCI,
           double r,
           double s):
    """get_ro(double Covar_ZvZv, np.ndarray [double, ndim=1] PCI, double r, double s)
    
    This function deduces the information effect coefficient, 
    conditional bias component, ro defined by the equations: 
    
        :math:`Covar (Zv,Zv^*) = \sum PCI^2 * r^n * s^n * ro^n`
        
    ro is necessary to account for the conditional bias in the 
    information effect.     
        
    The value of ro is deduced by finding the root of the equation 
    f_covar_ZvZv, using the classic Brent method (see 
    scipy.optimize.brentq)
    
    Parameters
    ----------
    r, s : float
        support effect and information effect (smoothing component)
    PCI : 1D array of floats
        hermite coefficients
    Covar_ZvZv : float
        Block covariance (correlation) between true Zv and estimate Zv* 
    
    Note
    ----
    :math:`Covar (Zv,Zv*) = var(Zv) - Kriging variance - LaGrange multiplier`
    
    see expression 7.47, page 97 on Basic Linear Geostatistics by 
    Margaret Armstrong.
    
    In case of information effect this can be calculated with a dummy 
    dataset in a single block representing future information, for 
    example blast holes. 
    
    Note that the slop of regression is  
    
    :math:`p = Covar (Zv,Zv^*) / (Covar (Zv,Zv^*)  - LaGrange multiplier)`
    
    """
    
    return brentq(f=f_covar_ZvZv, a=0, b=1, args=(s,r,PCI,Covar_ZvZv))


    
    
# ----------------------------------------------------------------------
#   Grade tonnage curves
# ----------------------------------------------------------------------
def gtcurve (cutoff, z, p, varred = 1, ivtyp = 0, zmin = None, zmax = None,
             ltail = 1, ltpar = 1, middle = 1, mpar = 1, utail = 1, utpar = 1,maxdis = 1000):
    """ gtcurve (cutoff, z, p, varred = 1, ivtyp = 0, zmin = None, zmax = None,
             ltail = 1, ltpar = 1, middle = 1, mpar = 1, utail = 1, utpar = 1,maxdis = 1000)
    
    Calculates grade and normalized tonnage (0 to 1) above a given cutoff 
    from input pairs of variable z, and probability p. The function calls the `pygslib.gslib.postik`
    to calculate probability and grade. 

    You can pass z and p corrected to block support with DGM, with postik parameter `varred == 1`,
    to get a DGM global change of support. Alternatively, you can pass z and p in point support
    and use a variance reduction factor  `varred < 1` and `ivol = 1` to do affine correction (`ivtyp=1`)
    or indirect lognormal correction (`ivtyp=2`). 

    Parameters
    ----------
    cutoff: 1D numeric array
        1D array with cutoff grades to plot
    z, p: 1D numeric array 
        variable (grade) and its probability (CDF). 
    varred: int, default  1
        variance reduction factor. If varred!=1 variance reduction will be 
        applied with affine correction (`ivtyp=1`), or lognormal correction (`ivtyp=2`)
    ivtyp: int, default  0
        type of change of support applied if varred!=1. 
        affine correction (`ivtyp=1`), and lognormal correction (`ivtyp=2`) are the two options. 
        for DGM global change of support preprocess z, p using `pygslib.nonlinear.anamor_blk`.
    zmin, zmax: floats, default None
        minimum and maximum of the z in the CDF. It can be outside the data interval. 
        If None will be set to input z minimum and maximum.
    ltail, ltpar, middle, mpar , utail, utpar: numeric
        lower, middle, and upper tail type and parameter
        ltail,mtail,utail=1 is linear interpolation.
        ltail,mtail,utail=2 is power model 
        ltail=3 is not implemented in this fuction, but availiable in `pygslib.gslib.postik`
    maxdis: integer, default 1000
        maximum discretization parameter

    Notes:
    ----
    see http://www.gslib.com/gslib_help/postik.html for postik gslib options

    """
    
    t = np.zeros([len(cutoff)])
    ga = np.zeros([len(cutoff)])
    gb = np.zeros([len(cutoff)])
    
    t[:] = np.nan
    ga[:] = np.nan
    gb[:] = np.nan
    
    ivol = 0
    
    if varred!=1.: 
        ivol = 1
        
    if zmin is None:
        zmin = min(z)
        
    if zmax is None:
        zmax = max(z)
    
    postik_parameters = {
            # output option, output parameter
            'iout'   : 2 ,   # 2 to calculate prob and mean above cutoff 
            'outpar' : None ,   # the cutoff. 
            # the thresholds
            'ccut1'  : z,   # these are Z values in the CDF
            # volume support?, type, varred
            'ivol'   : ivol,   # input int
            'ivtyp'  : ivtyp,   # input int
            'varred' : varred,   # input float
            # minimum and maximum Z value
            'zmin'   : zmin,   # input float
            'zmax'   : zmax,   # input float
            # lower,middle and upper tail: option, parameter
            'ltail'  : ltail,   # input int
            'ltpar'  : ltpar,   # input float
            'middle'  : middle,   # input int
            'mpar'  : mpar,   # input float
            'utail'  : utail,   # input int
            'utpar'  : utpar,   # input float
            # maximum discretization
            'maxdis' : maxdis,   # input int
            # 1D arrays with global distribution
            'vr'     : [0,0,0],   # input rank-1 array('f') with bounds (nc)
            'wt'     : [0,0,0],   # input rank-1 array('f') with bounds (nc)
            # 2D array with IK3D output (continuous)
            'p'      : [p]}   # input rank-2 array('f') with bounds (na,nc)
    
    for i in range(len(cutoff)):
        postik_parameters['outpar'] = cutoff[i]
        out1,out2,out3 = pygslib.gslib.postik(postik_parameters)
        
        t[i] = out1[0]
        ga[i] = out2[0]
        gb[i] = out3[0]
        
    return t,ga,gb    
    
def plotgt(cutoff, t, g, label, figsize = [6.4, 4.8]):
    """ plotgt(cutoff, t, g, label)
    
    Plots grade and tonnage above cutoff previously calculated

    Parameters
    ----------
    cutoff: 1D numeric array
        1D array with cutoff grades to plot
    t, g: 2D numeric array 
        this is an array of n,m array, where n is the number of tonnage curves
        and m = len(cutoff). 
    label: 1D array of strings
        the names of each grade and tonnage curve
    figsize: 1D array, default [6.4, 4.8]
        size of the matplotlib figue
    """
    
    # prepare ax and fig
    fig, ax1 = plt.subplots(figsize=figsize)
    ax1.set_xlabel('cutoff')
    ax1.set_ylabel('Tonnage')
    ax1.tick_params('y')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Grade')
    ax2.tick_params('y')
    
    # plot tonage
    for i,l in zip(t,label): 
        ax1.plot(cutoff, i, label = l)

    ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)     
        
    # plot grade
    for i in g:
        ax2.plot(cutoff, i)

    fig.tight_layout()
    plt.show()
    
    return fig    
    
    


# ----------------------------------------------------------------------
#   Uniform conditioning functions
# ----------------------------------------------------------------------

cpdef ucondit(float YV, 
              float yc,
              np.ndarray [double, ndim=1] PCI,
              float r=1., 
              float R=1., 
              float ro=1.): 
    
    """ ucondit(YV, yc, PCI, r=1., R=1., ro=1.)
    
    Computes uniform conditioning estimation in a panel with krided gaussian value YV,
    at cutoff yc, given support effect coefficients r, R, and information 
    effect coefficient ro

    Parameters
    ----------
    YV: numeric
        Kriged gaussian grade in a panel (get it from panel kriged values ZV using panel anamorphosis)
    PCI: 1D numeric array of double 
        PCI coefficients 
    yc: numeric
        cutoff in gaussian space (get it from zc using the SMU anamorphosis function)
    r,R,ro: numeric
        point and panel support effect coefficients, and information 
        effect coefficient
    
    Returns
    ----
    T, Q: 
        floats with tonnage and metal above cutoff

    Note
    ----
    Note that you may use gaussian panel grade (YV), and gaussian cutoff (yc) 
    as input. 
    
    You can get yc values from tabulated point transformations using the function 
    `pygslib.nonlinear.Z2Y_linear`. 

    YV can be interpolated directly in panels from y(x) on samples, or it can be 
    obtained transforming kriged panel values in non-gaussuan space (ZV) into 
    its gaussian equivalent using tabulated panel anamorphosis and the function 
    `pygslib.nonlinear.Z2Y_linear`

    """ 
    cdef float t, qn, qq
    cdef int K
    cdef float T 
    cdef float Q 
    cdef float M 
    cdef np.ndarray [double, ndim=1] HYV
    cdef np.ndarray [double, ndim=1] Hyc
    cdef np.ndarray [double, ndim=2] U
    
    # get general parameters 
    t=R/(r*ro)       # info and support effect (for no info make ro=1)
    K = PCI.shape[0]  # number of PCI
    
    # calculate tonnage
    T = 1- norm.cdf((yc-t*YV)/np.sqrt(1-t**2))  # this is ~ P[Zv>=zc] 
    

    # ERROR from here, the results are not OK

    # calculate metal (equation 3.17 from  C.T. Neufeld, 2015. Guide to Recoverable Reserves with Uniform Conditioning. Centre for Computational Geostatistics (CCG) Guidebook Series Vol. 4)
    # also C. Neufeld and C.V. Deutsch Calculating Recoverable Reserves with Uniform Conditioning
    HYV = recurrentH(np.array([YV]), K).ravel() # calculate hermite polynomials of YV
    Hyc = recurrentH(np.array([yc]), K).ravel() # calculate hermite polynomials of yc
    
    # get U using recurrence
    U = np.zeros([K,K])
    U[0,0] = 1-norm.cdf(yc)
    g = norm.pdf(yc)
    for k in range(1,K):
        U[0,k] = U[k,0] = -1/np.sqrt(k)*Hyc[k-1]*g

    for n in range(1,K):
        for p in range(n,K):
            U[n,p] = U[p,n]= -1/np.sqrt(n)*Hyc[p]*Hyc[n-1]*g + np.sqrt(p/n)*U[n-1][p-1]

    # Get Q
    Q = 0
    for n in range(0,K):
        qq = (t**n)*HYV[n]
        for p in range(0,K):
            qn = PCI[p]*(r**p)*(ro**p)*U[p][n]
            Q = Q + qq * qn                  # verify that Q = Q + PCI[p]*r**p  *(ro**p) *U[n][p]*t**n*H[n]; and enable info effect

    return T, Q
    

