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


# ----------------------------------------------------------------------
#   Transformation table
# ----------------------------------------------------------------------
cpdef ttable(z, w):
    """ttable(z, w)
    
    Creates a transformation table. 
    
    Parameters
    ---------
    z,w : 1D numpy arrays of floats 
        Variable and declustering weight.
    Returns:
    transin,transout : 1D numpy arrays with pairs of raw and gaussian values
    
    Note: 
    This function uses gslib.__dist_transf.ns_ttable
    """
    cdef int error
    
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
#   Transformation table
# ----------------------------------------------------------------------
cpdef stats(z, w, iwt = True, report = True):
    """stats(z, w)
    
    Reports some basic stats using declustering wights 
    
    Parameters
    ---------
    z,w : 1D numpy arrays of floats 
        Variable and declustering weight.
    iwt: boolean default True
        If True declustering weights will be used to calculate statsistics
    report : boolean default True
        If True a printout will be produced

    Returns:
    xmin,xmax, xcvr,xmen,xvar : floats
        minimum, maximum, coefficient of variation, media and variance
    
    Note: 
    This function uses gslib.__plot.probplt to optain the stats
    """
    cdef int error
    
    parameters_probplt = {
            'iwt'  : iwt,     #int, 1 use declustering weight
            'va'   : z,       # array('d') with bounds (nd)
            'wt'   : w}       # array('d') with bounds (nd), wight variable (obtained with declust?)



    binval,cl,xpt025,xlqt,xmed,xuqt,xpt975,xmin,xmax, xcvr,xmen,xvar,error = pygslib.gslib.__plot.probplt(**parameters_probplt)
    
    assert error < 1, 'There was an error = {} in the function gslib.__dist_transf.ns_ttable'.format(error)
    
    if report:
        print  'Stats Summary'
        print  '-------------'
        print  'Minimum        ', xmin
        print  'Maximum        ', xmax
        print  'CV             ', xcvr
        print  'Mean           ', xmen
        print  'Variance       ', xvar
        
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


cpdef ctrpoints(np.ndarray [double, ndim=1] z,
                np.ndarray [double, ndim=1] y,
                np.ndarray [double, ndim=1] ze):
    """ctrpoints(np.ndarray [double, ndim=1] z, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] ze)

    Infers some control points.
    
    The admissible minimum and maximum are just logic values defined 
    by user or by the nature of the problem, ej. grade values may be 
    between 0 and 100 percent. 

    The practical minimum and maximum are selected as follow

     - are within admissible interval
     - are within data minimum and maximum values
     - are within interval where ze do not fluctuates
    
    Parameters
    ----------
    z,y : 1D array of floats
        experimental transformation table 
    ze : 1D array of floats
        z values obtained with polynomial expansion of y, in other 
        ze may be calculated as expand_anamor(PCI,H,r=1) 

    Returns
    -------
    zamin,yamin,zamax,yamax,zpmin,ypmin,zpmax,ypmax : floats
        Control points
        
    """

    cdef float zamax, yamax, zamin, yamin
    cdef float zpmax, ypmax, zpmin, ypmin


    # initialize authorized interval as data max, min 
    zamax=z[y.shape[0]-1]
    yamax=y[y.shape[0]-1]
    zamin=z[0]
    yamin=y[0]

    # get the practical maximum interval
    zpmax=z[y.shape[0]-1]
    ypmax=y[y.shape[0]-1]
    for i in range (np.argmin(abs(y)), y.shape[0]-1):  # np.argmin(abs(y)) is around 50% prob. or y~0
        if ze[i]>ze[i+1] or ze[i]> zamax:
            zpmax=z[i-1]
            ypmax=y[i-1]
            break 

    # get the practical minimum interval
    zpmin=z[0]
    ypmin=y[0]
    for i in range(np.argmin(abs(y)), 1, -1):
        if ze[i-1]>ze[i] or ze[i]<zamin:
            zpmin=z[i+1]
            ypmin=y[i+1]
            break 

    return  zamin, yamin, zpmin, ypmin, zpmax, ypmax, zamax, yamax
     

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
    Z : 1D array of float64
        raw values corresponding to Y 
      
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
cpdef anamor(z, w, zmin, zmax, ltail=1, utail=1, ltpar=1, utpar=1, K=30):
    """anamor(z, w)
    """
    
    # a) get experimental transformation table
    transin, transout = ttable(z, w)
    
    # b) genarate a sequence of gaussian values
    gauss = np.linspace(-5,5,10000)
    
    # c) Get the back tansform using normal score transformation
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
    
    plt.plot(gauss,raw, '.')
    plt.plot(gauss,zana, '.')
    
    print 'Raw Variance', raw_var
    print 'Variance from PCI', PCI_var
    
    return PCI, raw, gauss, raw_var , PCI_var
    

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
    f_var_Zv, fit_PCI, scipy.optimize.brentq
    
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
#   Uniform conditioning functions
# ----------------------------------------------------------------------

#cpdef recurrentU(np.ndarray [double, ndim=2] H, float yc):
    
#    U =  H 
"""
cpdef ucondit(np.ndarray [double, ndim=1] ZV,
              np.ndarray [double, ndim=1] PCI, 
              float zc,
              float r=1., 
              float R=1., 
              float ro=1.): 
    
    
    # r block support, R panel support, ro info effect  
    
    
    cdef float t
    cdef int K
    cdef np.ndarray [double, ndim=1] T, Q, M
    cdef np.ndarray [double, ndim=2] H
    
    # get general parameters 
    t=R/(r*ro)       # info and support effect (for no info make ro=1)
    K = PCI.shape[0]
    yc = Z2Y_linear 
    YV = Z2Y_linear
    
    H = recurrentH(YV, K)
    
    T[:] = 1- norm.pdf((yc-t*YV)/np.sqrt(1-t**2))
    
    Q = np.zeros ([ZV.shape[0]])
    
    for i in range(K):
        for j in range(K): 
            Q = Q + t**i*H[i][:] * PCI[j]*r**j*ro**j
    
    #M = Q / T
    
    return T, Q, M
    
"""
