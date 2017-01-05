'''
PyGSLIB sandbox, Module to test new functions.  

Copyright (C) 2016 Adrian Martinez Vargas 

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


cimport numpy as np
import numpy as np
from scipy.interpolate import Rbf
import scipy.spatial.distance as scipy_dist
from scipy.interpolate import griddata




# ----------------------------------------------------------------------
#   Functions to interpolate smoth data
# ----------------------------------------------------------------------
def interpolator_UGdata(xt,yt,zt,x,y,z,v, method= 'linear', fill_value = .0):
    """
    Interpolates v from input 3D points (x,y,z) into target 
    locations (xt,yt,zt).  Points (x,y,z) may be unstructured data.
        	
	
    Parameters
    ----------
    xt,yt,zt : arrays
        target location
    x,y,z,v: arrays
        input coordinates and variable values
    method: str or callable, default('linear')
        one of {‘linear’, ‘nearest’}
    fill_value : float, default(.0)
        Value used to fill in for requested points outside of the 
        convex hull of the input points. 
    
    
    Returns 
    -------
    vt : array floats
        values interpolated at locations xt,yt,zt
        
    Notes
    -------
    Using scipy.interpolate.griddata 
    
    The linear method is the barycentric interpolation on triangulated input data with Qhull. 
    
    """
    
    # check if there are duplicates
    # assert scipy_dist.pdist(np.array([x,y,z]).T).min()  > mindis, 'Error: minimum distance in x,y,y <= {}'.format(mindis)
          
    vt =  griddata(np.array([x,y,z]).T, v, np.array([xt,yt,zt]).T, method=method, fill_value = fill_value)
    
    
    return  vt



def interpolator_RBF3D(xt,yt,zt,x,y,z,v, function= 'linear', smooth = .0, mindis = 0.01):
    """
    Interpolates v from input points (x,y,z) into target 
    locations (xt,yt,zt) using Radial Basis Functions (RBF)
    
    
    Parameters
    ----------
    xt,yt,zt : arrays
        target location
    x,y,z,v: arrays
        input coordinates and variable values
    function: str or callable, default('multiquadric')
        'multiquadric': sqrt((r/self.epsilon)**2 + 1)
        'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
        'gaussian': exp(-(r/self.epsilon)**2)
        'linear': r
        'cubic': r**3
        'quintic': r**5
        'thin_plate': r**2 * log(r)
    smooth: float, default(.0)
        Values greater than zero increase the smoothness of the 
        approximation. 0 is for interpolation (default), the function 
        will always go through the nodal points in this case.
    
    
    Returns 
    -------
    vt : array floats
        values interpolated at locations xt,yt,zt
        
    Notes
    -------
    Using RBF from scipy. We recommend to see RBF documentation
    at https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
    
    
    TODO
    ------
    Move this data to a new file / module
    
    Example
    -------
    
    >>> vt = pygslib.sandbox.interpolator_RBF3D(xt,yt,zt,x,y,z,v, function= 'linear', smooth = .0)
    
    
    """
    
    # check if there are duplicates
    assert scipy_dist.pdist(np.array([x,y,z]).T).min()  > mindis, 'Error: minimum distance in x,y,y <= {}'.format(mindis)
    
    
    rbfi = Rbf(x, y, z, v, function = function, smooth =smooth)  # radial basis function interpolator instance
    
    vt = rbfi(xt, yt, zt)   # interpolated values vt
    
    return  vt

def interpolator_RBF2D(xt,yt,x,y,v, function= 'linear', smooth = .0, mindis = 0.01):
    """
    Interpolates v from input points (x,y) into target 
    locations (xt,yt) using Radial Basis Functions (RBF)
    
    
    Parameters
    ----------
    xt,yt : arrays
        target location
    x,y,v: arrays
        input coordinates and variable values
    function: str or callable, default('multiquadric')
        'multiquadric': sqrt((r/self.epsilon)**2 + 1)
        'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
        'gaussian': exp(-(r/self.epsilon)**2)
        'linear': r
        'cubic': r**3
        'quintic': r**5
        'thin_plate': r**2 * log(r)
    smooth: float, default(.0)
        Values greater than zero increase the smoothness of the 
        approximation. 0 is for interpolation (default), the function 
        will always go through the nodal points in this case.
    
    
    Returns 
    -------
    vt : array floats
        values interpolated at locations xt,yt
        
    Notes
    -------
    Using RBF from scipy. We recommend to see RBF documentation
    at https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
    
    
    TODO
    ------
    Move this data to a new file / module
    
    Example
    -------
    
    >>> vt = pygslib.sandbox.interpolator_RBF2D(xt,yt,x,y,v, function= 'linear', smooth = .0)
    
    
    """
    
    # check if there are duplicates
    assert scipy_dist.pdist(np.array([x,y]).T).min()  > mindis, 'Error: minimum distance in x,y <= {}'.format(mindis)
    
    
    
    rbfi = Rbf(x, y, v, function = function, smooth =smooth)  # radial basis function interpolator instance
    
    vt = rbfi(xt, yt)   # interpolated values vt
    
    return  vt

# ----------------------------------------------------------------------
#   Functions to warp/unfold data
# ----------------------------------------------------------------------

def Unfold3D(x1,y1,z1,x0,y0,z0,xd0,yd0,zd0,function= 'linear', smooth = .0, interpolator = 'UGdata', fill_value = np.nan):
    """
    
    This function unfolds data at coordinates Xd using a set of control 
    points/lines with folded coordinates X0 and unfolded coordinates X1.  
    
    To unfold a warp vector is calculated at control point as 
    ``W = X1-X0``, then the ``W`` vectors are interpolated to each 
    location ``Xd0`` using the function ``interpolator_RDB()`` to obtain
    ``Wd``. The unfolded coordinates are calculated as ``Xd1 = Xd0 + Wd``
        
    
    Parameters
    ----------
    x0,y0,z0 : arrays
        control points at folded space
    x1,y1,z1 : arrays
        control points at unfolded space
    xd,yd,zd : arrays
        target location to be unfolded
    function: str or callable, default('multiquadric')
        'multiquadric': sqrt((r/self.epsilon)**2 + 1)
        'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
        'gaussian': exp(-(r/self.epsilon)**2)
        'linear': r
        'cubic': r**3
        'quintic': r**5
        'thin_plate': r**2 * log(r)
		'nearest' : (only valid if using interpolator 'UGdata')
    smooth: float, default(.0)
        Values greater than zero increase the smoothness of the 
        approximation. 0 is for interpolation (default), the function 
        will always go through the nodal points in this case.
    interpolator : str, default('UGdata')
		one of {'UGdata', 'RBF3D'}. This is the interpolator used to 
		interpolate desplacement vectors wx, wy, wz. 
	fill_value : float, default(np.nan)
		if using interpolator == 'UGdata' fill_value is the value 
		asigned outside of the convex hull of the input points.
    
    Returns 
    -------
    xd1,yd1,zd1 : arrays of floats
        coordinates of data points in unfolded space
        
    wx, wy, wz : arrays of floats
        desplacement vectors
        
    Notes
    -------
    This functions uses RBF from scipy. We recommend to see RBF documentation
    at https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
    
    Unfold3D interpolates with 3D RBF or linear barycentric interpolation
	
    
    TODO
    ------
    Move this data to a new file / module
    
    """

    # calculate warp vector at control points
    wx=x1-x0
    wy=y1-y0
    wz=z1-z0
    
    
    # interpolate warp vectors at data points
    if interpolator == 'RBF3D': 
    
        # we use RBF 

        wdx = interpolator_RBF3D(xd0, yd0, zd0, x0, y0, z0, wx,
                            function= function,
                            smooth = smooth)
        wdy = interpolator_RBF3D(xd0, yd0, zd0, x0, y0, z0, wy,
                            function= function,
                            smooth = smooth)
        wdz = interpolator_RBF3D(xd0, yd0, zd0, x0, y0, z0, wz,
                            function= function,
                            smooth = smooth)
    elif interpolator == 'UGdata':                   
                        
        wdx = interpolator_UGdata(xd0, yd0, zd0, x0, y0, z0, wx, 
                               method= function, fill_value = fill_value)

        wdy = interpolator_UGdata(xd0, yd0, zd0, x0, y0, z0, wy, 
                               method= function, fill_value = fill_value)

        wdz= interpolator_UGdata(xd0, yd0, zd0, x0, y0, z0, wz,  
                               method= function, fill_value = fill_value)

    # unfold
    xd1 = xd0 + wdx
    yd1 = yd0 + wdy
    zd1 = zd0 + wdz
    
    return xd1, yd1, zd1, wx, wy, wz
    
    
# ----------------------------------------------------------------------
#   Functions to warp/unfold data
# ----------------------------------------------------------------------
def Unfold2D(x1,y1,z1,x0,y0,z0,xd0,yd0,zd0,function= 'linear', smooth = .0):
    """
    
    This function unfolds data at coordinates Xd using a set of control 
    points/lines with folded coordinates X0 and unfolded coordinates X1.  
    
    To unfold a warp vector is calculated at control point as 
    ``W = X1-X0``, then the ``W`` vectors are interpolated to each 
    location ``Xd0`` using the function ``interpolator_RDB()`` to obtain
    ``Wd``. The unfolded coordinates are calculated as ``Xd1 = Xd0 + Wd``
        
    
    Parameters
    ----------
    x0,y0,z0 : arrays
        control points at folded space
    x1,y1,z1 : arrays
        control points at unfolded space
    xd,yd,zd : arrays
        target location to be unfolded
    function: str or callable, default('multiquadric')
        'multiquadric': sqrt((r/self.epsilon)**2 + 1)
        'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
        'gaussian': exp(-(r/self.epsilon)**2)
        'linear': r
        'cubic': r**3
        'quintic': r**5
        'thin_plate': r**2 * log(r)
    smooth: float, default(.0)
        Values greater than zero increase the smoothness of the 
        approximation. 0 is for interpolation (default), the function 
        will always go through the nodal points in this case.
    
    
    Returns 
    -------
    xd1,yd1,zd1 : arrays
        coordinates of data points in unfolded space

    wx, wy, wz : arrays of floats
        desplacement vectors

        
    Notes
    -------
    This functions uses RBF from scipy. We recommend to see RBF documentation
    at https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
    
    Unfold2D interpolates with 2D RBF in plane XY
    
    TODO
    ------
    Move this data to a new file / module
    
    """

    # calculate warp vector at control points
    wx=x1-x0
    wy=y1-y0
    wz=z1-z0
    
    
    # interpolate warp vectors at data points

    wdx = interpolator_RBF2D(xd0, yd0,  x0, y0, wx,
                        function= function,
                        smooth = smooth)
    wdy = interpolator_RBF2D(xd0, yd0, x0, y0, wy,
                        function= function,
                        smooth = smooth)
    wdz = interpolator_RBF2D(xd0, yd0, x0, y0, wz,
                        function= function,
                        smooth = smooth)

    # unfold
    xd1 = xd0 + wdx
    yd1 = yd0 + wdy
    zd1 = zd0 + wdz
    
    return xd1, yd1, zd1, wx, wy, wz
