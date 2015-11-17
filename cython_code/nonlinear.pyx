'''
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
'''

'''
Code based on paper:

A Step by Step Guide to Bi-Gaussian Disjunctive Kriging, by 
Julian M. Ortiz, Bora Oz, Clayton V. Deutsch
Geostatistics Banff 2004
Volume 14 of the series Quantitative Geology and Geostatistics pp 1097-1102

see also:

http://www.ccgalberta.com/ccgresources/report05/2003-107-dk.pdf
http://www.ccgalberta.com/ccgresources/report04/2002-106-hermite.pdf
http://www.ccgalberta.com/ccgresources/report06/2004-112-inference_under_mg.pdf
Mining Geostatistics: A. G. Journel, Andre G. Journel, C. J
'''


cimport numpy as np
import numpy as np


# ----------------------------------------------------------------------
#   Functions for punctual gaussian anamorphosis 
# ----------------------------------------------------------------------

#the recurrent formula for normalized polynomials
cpdef recurrentH(np.ndarray [double, ndim=1] y, K=30):
    """ 
    recurrentH(y, K=30)
    
    Calculate the hermite polynomials with the recurrent formula
    
    Parameters
    ----------
    y : 1D array of floats
        Gaussian values calculated for the right part of the bin.
    K  : integer, default 30
        Number of hermite polynomials 

    Returns
    -------
    H : 2D array of floats
        Hermite monomials H(i,y) with shape [K+1,len(y)]
      
    See Also
    --------
    pygslib.__dist_transf.anatbl
       
    Notes
    -----  
    The y values may be calculated on the right side of the bin, 
    as shown in fig VI.13, page 478 of Mining Geostatistics: 
    A. G. Journel, Andre G. Journel, C. J. The function 
    pygslib.__dist_transf.anatbl was prepared to provide these values,
    considering declustering weight if necessary. 
    
    The results from pygslib.__dist_transf.ns_ttable are inappropriate 
    for this calculation because are the mid interval in the bin.  
    """
    assert(K>=1)
    
    cdef np.ndarray [double, ndim=2] H
    cdef int k
    
    H=np.ones((K+1,len(y))) 
    H[0,:]=1                 #first monomial
    H[1,:]=-y               #second monomial
    
    # recurrent formula
    for k in range(1,K):
        H[k+1,:]= -1/np.sqrt(k+1)*y*H[k,:]-np.sqrt(k/float(k+1))*H[k-1,:]
    
    return H   #this is a 2D array of H (ki,yi)


#fit PCI for f(x)=Z(x)
cpdef fit_PCI(z,y,H):
    """ 
    fit_PCI(z,y,H)
    
    Fit the hermite coefficient (PCI) 
    
    Parameters
    ----------
    z  : 1D array of floats
        Raw values sorted
    y  : 1D array of floats
        Gaussian values calculated for the right part of the bin.
    K  : integer, default 30
        Number of hermite polynomials 

    Returns
    -------
    H : 2D array of floats
        Hermite monomials H(i,y) with shape [K+1,len(y)]
      
    See Also
    --------
    recurrentH
       
    Notes
    -----  
    
    """
    PCI=np.zeros(H.shape[0])
    PCI[0]=np.mean(z)
    
    for p in range(1,H.shape[0]):
        for i in range(1,H.shape[1]):
            PCI[p]+=(z[i-1]-z[i])*1/np.sqrt(p)*H[p-1,i]*norm.pdf(y[i])
    
    return PCI



#get variance from PCI
cpdef var_PCI(PCI):
    """
    put comment here
    """
    a=PCI[1:]**2
    return np.sum(a)

#expand anamorphosis
cpdef expand_anamor(PCI,H,r=1,ro=1):
    """
    from now ising the same data y
    to do use arbitrary data y
    
    If there is no conditional bias, then s=r*ro. 
    This is equivalent to say that the blocks are equivalent 
    to larger blocks with perfect estimation "perfect know grade". 

    """
    Z=np.zeros(H.shape[1])
    for i in range(len(PCI)):
        Z+=PCI[i]*H[i,:]*(r*ro)**i
    
    return Z
