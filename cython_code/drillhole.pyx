'''
PyGSLIB Drillhole, Module to handle drillhole data, desurvey 
interval tables and other drillhole relate process.  

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

import vtk
cimport numpy as np
import numpy as np
from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport acos
from libc.math cimport asin
from libc.math cimport tan
from libc.math cimport atan2
import pandas as pd
import warnings

#-------------------------------------------------------------------
#  General functions for desurvey 
#-------------------------------------------------------------------
cpdef ang2cart( float azm,
                float dip):
    """
    ang2cart(azm, dip)
    
    Convert azimuth and dip to x, y, z. 
    
    Return the x,y,z coordinates of a 1 unit vector with origin of 
    coordinates at 0,0,0 and direction defined by azm (azimuth) and 
    dip (downward positive) angles. 
    
    
    Parameters
    ----------
    azm, dip : float, in degrees
        The direction angles azimuth, with 0 or 360 pointing north and
        the dip angle measured from horizontal surface positive downward
    
    Returns
    -------
    out : tuple of floats, ``(x, y, z)``
        Cartesian coordinates.
    
    See Also
    --------
    cart2ang
    
    Notes
    -----
    This function is to convert direction angle into a cartesian 
    representation, which are easy to interpolate. To convert
    back x,y,z values to direction angles use the function cart2ang.
    
    Examples
    --------
    
    >>> ang2cart(45,75)
    (0.1830127239227295, 0.1830127090215683, -0.9659258127212524)

    """

    # output
    cdef: 
        float x
        float y
        float z

    # internal
    cdef:
        float razm
        float rdip
        float DEG2RAD

    DEG2RAD=3.141592654/180.0

    # convert degree to rad and correct sign of dip
    razm = azm * DEG2RAD
    rdip = -dip * DEG2RAD

    # do the conversion
    x = sin(razm) * cos(rdip)
    y = cos(razm) * cos(rdip)
    z = sin(rdip)

    return x,y,z



cpdef cart2ang( float x,
                float y,
                float z):
    """
    cart2ang(x, y, z)
    
    Convert x, y, z to azimuth and dip. 
    
    Return the azimuth and dip of a 1 unit vector with origin of 
    coordinates at p1 [0,0,0] and p2 [x, y, z]. 
    
    
    Parameters
    ----------
    x, y, z : float
        Coordinates x, y, z of one unit length vector with origin of 
        coordinates at p1 [0,0,0]
        
    Warning
    -------
    If x, y or z are outside the interval [-1,1] the values will be 
    truncated to 1. If the coordinates are outside this 
    interval due to rounding error then result will be ok, otherwise
    the result will be wrong 
    
    Returns
    -------
    out : tuple of floats, ``(azm, dip)``
        The direction angles azimuth, with 0 or 360 pointing north and
        the dip angle measured from horizontal surface positive downward
    
    See Also
    --------
    ang2cart
    
    Notes
    -----
    See function cart2ang.
    
    Examples
    --------
    
    >>> cart2ang(0.18301, 0.18301, -0.96593)
    (45.0, 75.00092315673828)
    
    """

    # out 
    cdef:
        float azm
        float dip


    # internal
    cdef:
        float razm
        float rdip
        float RAD2DEG
        float pi
        float EPSLON=1.0e-4

    if abs(x)+EPSLON>1.001 or abs(y)+EPSLON>1.001 or abs(z)+EPSLON>1.001:
        #print x,y,z
        warnings.warn('cart2ang: a coordinate x, y or z is outside the interval [-1,1]')
        
    # check the values are in interval [-1,1] and truncate if necessary
    if x>1.: x=1.
    if x<-1.: x=-1.
    if y>1.: y=1.
    if y<-1.: y=-1.   
    # in case z>1 or z<-1 asin will be out of domain and will produce an error 
    if z>1.: z=1.
    if z<-1.: z=-1.

    RAD2DEG=180.0/3.141592654
    pi = 3.141592654

    if x!=0. and y!= 0.: 
        azm= atan2(x,y)
        if azm<0.:  
            azm= azm + pi*2
        azm = azm * RAD2DEG
    else: 
        azm = 0.


    dip = -asin(z) * RAD2DEG
     
    return azm, dip


cpdef interp_ang1D( float azm1,
                    float dip1,
                    float azm2,
                    float dip2,
                    float len12,
                    float d1):
    """    
    interp_ang1D(azm1, dip1, azm2, dip2, len12, d1)
    
    Interpolate the azimuth and dip angle over a line. 
    
    Given a line with length ``len12`` and endpoints with direction 
    angles ``azm1, dip1, azm2, dip2``, this function returns the 
    direction angles of a point over this line, located at a distance 
    d1 of the point with direction angles ``azm1, dip1``. 
    
    
    Parameters
    ----------
    azm1, dip1, azm2, dip2, len12, d1 : float
        azm1, dip1, azm2, dip2 are direction angles azimuth, with 0 or 
        360 pointing north and dip angles measured from horizontal 
        surface positive downward. All these angles are in degrees.
        len12 is the length between a point 1 and a point 2 and 
        d1 is the positive distance from point 1 to a point located 
        between point 1 and point 2. 
        
        
    Returns
    -------
    out : tuple of floats, ``(azm, dip)``
        The direction angles azimuth, with 0 or 360 pointing north and
        the dip angle measured from horizontal surface positive downward
    
    See Also
    --------
    ang2cart, cart2ang
    
    Notes
    -----
    The output direction angles are interpolated using average weighted  
    by the distances ``d1`` and ``len12-d1``. To avoid issues with 
    angles the azimuth and dip are converted to x,y,z, then interpolated 
    and finally converted back to azimuth,dip angles.
    
    Example
    --------
    
    >>> interp_ang1D(azm1=45, dip1=75, azm2=90, dip2=20, len12=10, d1=5)
    (80.74163055419922, 40.84182357788086)

    """
    # output
    cdef:
        float azm
        float dip

    # internal
    cdef:
        float x1
        float y1
        float z1
        float x2
        float y2
        float z2
        float x
        float y
        float z


    # convert angles to coordinates
    x1,y1,z1 = ang2cart(azm1,dip1)
    x2,y2,z2 = ang2cart(azm2,dip2)

    # interpolate x,y,z
    x = x2*d1/len12 + x1*(len12-d1)/len12 
    y = y2*d1/len12 + y1*(len12-d1)/len12
    z = z2*d1/len12 + z1*(len12-d1)/len12

    # get back the results as angles
    azm,dip = cart2ang(x,y,z)

    return azm, dip

cpdef dsmincurb( float len12,
                 float azm1,
                 float dip1,
                 float azm2,
                 float dip2):
                     
    """    
    dsmincurb(len12, azm1, dip1, azm2, dip2)
    
    Desurvey one interval with minimum curvature 
    
    Given a line with length ``len12`` and endpoints p1,p2 with 
    direction angles ``azm1, dip1, azm2, dip2``, this function returns 
    the differences in coordinate ``dz,dn,de`` of p2, assuming
    p1 with coordinates (0,0,0)
    
    Parameters
    ----------
    len12, azm1, dip1, azm2, dip2: float
        len12 is the length between a point 1 and a point 2.
        azm1, dip1, azm2, dip2 are direction angles azimuth, with 0 or 
        360 pointing north and dip angles measured from horizontal 
        surface positive downward. All these angles are in degrees.
        
        
    Returns
    -------
    out : tuple of floats, ``(dz,dn,de)``
        Differences in elevation, north coordinate (or y) and 
        east coordinate (or x) in an Euclidean coordinate system. 
    
    See Also
    --------
    ang2cart, 
    
    Notes
    -----
    The equations were derived from the paper: 
        http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
    
    The minimum curvature is a weighted mean based on the
    dog-leg (dl) value and a Ratio Factor (rf = 2*tan(dl/2)/dl )
    if dl is zero we assign rf = 1, which is equivalent to  balanced 
    tangential desurvey method. The dog-leg is zero if the direction 
    angles at the endpoints of the desurvey intervals are equal.  
    
    Example
    --------
    
    >>> dsmincurb(len12=10, azm1=45, dip1=75, azm2=90, dip2=20)
    (7.207193374633789, 1.0084573030471802, 6.186459064483643)

    """

    # output
    cdef:
        float dz
        float dn
        float de 


    # internal 
    cdef:
        float i1
        float a1
        float i2
        float a2
        float DEG2RAD
        float rf
        float dl

    DEG2RAD=3.141592654/180.0

    i1 = (90 - dip1) * DEG2RAD
    a1 = azm1 * DEG2RAD

    i2 = (90 - dip2) * DEG2RAD
    a2 = azm2 * DEG2RAD

    # calculate the dog-leg (dl) and the Ratio Factor (rf)
    dl = acos(cos(i2-i1)-sin(i1)*sin(i2)*(1-cos(a2-a1))) 

    if dl!=0.: 
        rf = 2*tan(dl/2)/dl  # minimum curvature
    else:
        rf=1                 # balanced tangential



    dz = 0.5*len12*(cos(i1)+cos(i2))*rf
    dn = 0.5*len12*(sin(i1)*cos(a1)+sin(i2)*cos(a2))*rf
    de = 0.5*len12*(sin(i1)*sin(a1)+sin(i2)*sin(a2))*rf

    return dz,dn,de

cpdef desurv1dh(int indbs,
               int indes,
               np.ndarray[double, ndim=1] ats,
               np.ndarray[double, ndim=1] azs,
               np.ndarray[double, ndim=1] dips,
               float xc,
               float yc,
               float zc,
               float lpt,
               bint warns=True):

    """
    desurv1dh(indbs, indes, ats, azs, dips, xc, xy, zc, lpt)
    
    Desurvey one point with minimum curvature given survey array
    and collar coordinates
    
    It takes an array of survey points ``ats, azs, dips`` with a given
    drillhole starting and ending at index ``indbs,indes``, the 
    coordinates of the collar xc, xy, zc and derurvey a point at depth 
    ``lpt`` (any point at interval table, ex. a FROM interval in assay). 
    
    Warning
    --------
    The survey may have ats==0 at index indbs but this is not 
    checked at this function. 
    
    
    Parameters
    ----------
    len12, azm1, dip1, azm2, dip2: float
        len12 is the length between a point 1 and a point 2.
        azm1, dip1, azm2, dip2 are direction angles azimuth, with 0 or 
        360 pointing north and dip angles measured from horizontal 
        surface positive downward. All these angles are in degrees.
        
        
    Returns
    -------
    out : tuple of floats, ``(azt,dipt,xt,yt,zt)``
        Direction angles and coordinates at a point lpt
    
    See Also
    --------
    Drillhole.desurvey 
    
    Notes
    -----
    This function is a convenience function call by Drillhole.desurvey
    it was exposed here for validation purpose. 
    
    TODO
    ----
    Make this function hidden. 
    
   
    """
    
    # output (angles at begin, mid and end interval)
    cdef:
        float azt
        float dipt
        float xt
        float yt
        float zt

    # internal
    cdef:
        int i
        int j
        float a,
        float b
        float azm1
        float dip1
        float azm2
        float dip2
        float len12
        float d1
        float EPSLON=1.0e-4
        float xa
        float ya
        float za
        float xb
        float yb
        float zb
        float dz
        float dn
        float de

    assert ats[indbs]<EPSLON, 'first survey > 0 at %d' % indbs

    for i in range (indbs,indes):
        
        # get the segment [a-b] to test interval
        a=ats[i]
        b=ats[i+1]
        azm1 = azs[i]
        dip1 = dips[i]
        azm2 = azs[i+1]
        dip2 = dips[i+1]
        len12 = ats[i+1]-ats[i]
        # desurvey at survey table
        if ats[i]>=-EPSLON and ats[i]<=EPSLON: #zero interval?
            xa = xc
            ya = yc
            za = zc
            # desurvey interval at survey... table
            
            dz,dn,de = dsmincurb(len12,azm1,dip1,azm2,dip2)

            xb = xa+de
            yb = ya+dn
            zb = za-dz
            
        else:
            xa = xb
            ya = yb
            za = zb
            # desurvey interval at zurvey... table
            
            dz,dn,de = dsmincurb(len12,azm1,dip1,azm2,dip2)
            xb = xa+de
            yb = ya+dn
            zb = za-dz


        # test if we are in the interval, interpolate angles
        if lpt>=a  and lpt<b:
            d1= lpt- a

            azt,dipt = interp_ang1D(azm1,dip1,azm2,dip2,len12,d1)

            # desurvey at interval table 
            dz,dn,de = dsmincurb(d1,azm1,dip1,azt,dipt)

            xt = xa+de
            yt = ya+dn
            zt = za-dz 

            return azt, dipt, xt, yt, zt
    
    a=ats[indes]
    azm1 = azs[indes]
    dip1 = dips[indes]
    azt = azm1
    dipt = dip1
    # the point is beyond the last survey? 
    if lpt>=a:
        d1= lpt- a
        # desurvey at interval table 
        dz,dn,de = dsmincurb(d1,azm1,dip1,azt,dipt)
        xt = xb+de
        yt = yb+dn
        zt = zb-dz 
        
        if warns==True: 
            warnings.warn('\n point beyond the last survey point at %s' % indes)
    else:
        if warns==True: 
            warnings.warn('\n not interval found at survey, at %s' % indes)

    return   azt, dipt, xt, yt, zt


#-------------------------------------------------------------------
#  General functions for compositing 
#-------------------------------------------------------------------
cpdef composite1dh(double[:] ifrom, 
                   double[:] ito, 
                   double[:] ivar, 
                   double cint = 1., 
                   double minlen=-1.):
    '''
    cfrom, cto, clen, cvar, cacum= composite(ifrom, ito, ivar, cint = 1, minlen=-1)
    
    Composite intervals in a single drillhole. The From-To intervals may be sorted. 
    
    Note
    ----
    See explanation of this implementation at 
    http://opengeostat.com/downhole-compositing-algorithm/
       
    
    Parameters
    ----------
    ifrom, ito:     1D arrays of floats 
        From - To  interval
    ivar :   1D array of floats
        variable to be composited
    cint: Optional, float, default 1.  
        length of the compositing intervals
    minlen: Optional, float, default -1.  
        minimum length of the composite, if <=0 then minlen = cint/2.
    
    Return
    ------- 
    (cfrom, cto, clen, cvar, cacum)
    cfrom, cto:  1D arrays of floats
         From, To composited intervals 
    clen, cvar, cacum:  1D arrays of floats     
         total length of intervals composited
         variable composited 
         variable accumulated
    
    '''
    
    
    assert len(ifrom) == len(ito)==len(ivar), "Error: ifrom, ito or ivar with different shapes"
    assert cint>0, "Error compositing length (cint) <= 0 "
    #assert all(ifrom < ito), "Error: ifrom >= ito, wrong or zero length intervals"
    #assert all(np.isfinite(ifrom)),"Error: ifrom with not finite elements"
    #assert all(np.isfinite(ito)),"Error: ito with not finite elements"
    
    cdef long i, l, k, ncomp, nintrb
    

    #get some array parameters     
    if minlen<0:
        minlen= cint/2.

    ncomp = int(ito[-1]/cint + 1)  # number of composites 
    nintrb = len(ifrom)           # number of intervals


    # define output arrays
    np_cfrom = np.zeros(ncomp)     # from 
    np_cto = np.zeros(ncomp)       # to
    np_clen = np.zeros(ncomp)      # length of composited information, may be != cto - cfrom
    np_cvar = np.zeros(ncomp)      # composited variable
    np_cvar[:] = np.nan
    np_cacum = np.zeros(ncomp)     # accumulate variable in the interval (useful for category) 
    np_iprop =np.zeros(nintrb)     # size of interval within composite
    np_icum =np.zeros(nintrb)      # size of interval within composite


    #declare the composite/internal arrays (memory slice)
    cdef:
        double[::1] cfrom = np_cfrom
        double[::1] cto = np_cto
        double[::1] clen = np_clen
        double[::1] cvar = np_cvar
        double[::1] cacum = np_cacum
        double[::1] iprop = np_iprop
        double[::1] icum = np_icum

    #update some from to intervals
    #TODO: at this point add constraints for contacts
    for i in range(ncomp):
        cfrom[i] = i*cint
        cto[i] =  (i+1)*cint

    # for each composite 
    for i in range (ncomp): 

        # initialize proportions 
        iprop[:]=0
        icum[:]=0
        
        #for each interval 
        for l in range (nintrb):
            
            # ignore interval if variable is nan 
            if np.isnan(ivar[l]):
                continue

            # case a, below the composite 
            if ifrom[l]>=cto[i]: 
                break

            # case b, over the composite
            if ito[l]<=cfrom[i]: 
                continue

            # --these are overlap--

            # case overlap top or contained
            if ito[l]>cfrom[i] and ito[l]<=cto[i]:     

                # case c, interval in composite
                if ifrom[l]>=cfrom[i]: 
                    iprop[l] = ito[l]-ifrom[l]
                    icum[l] = ivar[l]*iprop[l]


                # case d, overlap top
                else:
                    iprop[l] = ito[l]-cfrom[i]
                    icum[l] = ivar[l]*iprop[l]


            # case e, composite in interval
            if ifrom[l]<cfrom[i] and ito[l]>cto[i]: 
                iprop[l] = cto[i]-cfrom[i]
                icum[l] = ivar[l]*iprop[l]
                continue


            # case f, overlap bottom
            if ifrom[l]>=cfrom[i] and ifrom[l]<cto[i] and ito[l]>cto[i]:
                iprop[l] = cto[i]-ifrom[l]
                icum[l] = ivar[l]*iprop[l]
                continue

        clen[i] = np.nansum(iprop)
        
        
        if clen[i]>= minlen:  
            cacum[i] = np.nansum(icum)
            cvar[i] = cacum[i]/clen[i]  # wighted average
        #else: 
        #    cvar[i] = np.nan
        #    cacum[i] = np.nan

    # this algorithm return an empty interval if there is no residual
    # to fix this we remove the last interval if len < minlen
    # this also apply the minlen constraint
            
    return np_cfrom[np_clen>= minlen], np_cto[np_clen>= minlen], np_clen[np_clen>= minlen], np_cvar[np_clen>= minlen], np_cacum[np_clen>= minlen]

#-------------------------------------------------------------------
#  General functions for fill gaps and merge
#-------------------------------------------------------------------
cdef min_int(double la, 
             double lb, 
             double ia, 
             double ib, 
             double tol=0.01):
    """
    ia, ib, l = min_int(la, lb, ia, ib, tol=0.01)
    
    This is a function to be used from merge_one_dhole
    
    Given two complete drillholes A, B (no gaps and up to the end of the drillhole), this function returns
    the smaller of two intervals FromA[ia] FromB[ib] and updates the 
    indices ia and ib. There are three posible outcomes
    
    - FromA[ia]==FromB[ib]+/- tol. Returns mean of FromA[ia], FromB[ib] and ia+1, ib+1
    - FromA[ia] <FromB[ib]. Returns FromA[ia] and ia+1, ib
    - FromA[ia]> FromB[ib]. Returns FromB[ia] and ia, ib+1
    
    
    
    """
    
    # equal ?
    if (lb-1)<=la<=(lb+1):
        return ia+1, ib+1, (la+lb)/2
        
    # la < lb ?
    if la<lb:
        return ia+1, ib, la
    
    # lb < la ?
    if lb<la:
        return ia, ib+1, lb

cdef merge_one_dhole(double[:] la,
              double[:] lb, 
              long[:] ida, 
              long[:] idb, 
              double tol=0.01):
    """
    n, lab, newida, newidb = merge_one_dhole(la,lb, ida, idb, tol=0.01)
    
    Function to merge one drillhole. 
    
    Notes: 
    - Full drillhole is required (no gaps) and note that length are used 
    instead From-To. 
    - The length at the end of the drillhole must be the same
    
    Returns:
    n: number of intervals
    lab: length from collar of each interval
    newida, newidb: Id linking data with the original tables a, b.
    
    """
    
    # general declarations
    cdef:
        int ia=0
        int ib=0
    
        int maxia= len (la)
        int maxib= len (lb)
        int maxiab= len (lb) + len (la)
        bint inhole = True
        int n=-1 
    
    # prepare output as numpy arrays
    np_newida= np.zeros(maxiab, dtype=int)
    np_newidb= np.zeros(maxiab, dtype=int)
    np_lab= np.zeros(maxiab, dtype=float)
    
    # get memory view of the numpy arrays    
    cdef:
        long[::1] newida = np_newida
        long[::1] newidb = np_newidb
        double[::1] lab = np_lab
    
    #check those are complete dholes
    assert la[maxia-1]==lb[maxib-1]
    
    #loop on drillhole
    while inhole: 
        # get the next l interval and l idex for drillhole a and b
        ia, ib, l = min_int(la[ia], lb[ib], ia, ib, tol=0.01)
        n+=1
        newida[n]=ida[ia-1]
        newidb[n]=idb[ib-1]
        lab[n]=l
        
        #this is the end of hole (this fails if maxdepth are not equal)
        if ia==maxia or ib==maxib: 
            inhole=False
    
    return n, np_lab[:n+1], np_newida[:n+1], np_newidb[:n+1]

cdef fillgap1Dhole(double[:] in_f, 
            double[:] in_t, 
            long[:] id, 
            double tol=0.01,
            double endhole=-1):
    """
    np_nf,np_nt,np_nID,np_gap,np_overlap = illgap(in_f, in_t, id, tol=0.01, endhole=-1)
    
    Function to fill gaps in one drillhole. 
    
    Input:
    in_f, in_t, id: from, to intervals and interval ID.
        The insterval ID is required to link back sample values after
        adding gaps
    tol: default 0.01. Tolerance
        gaps and overlaps <= tolerance will be ignored
    endhole: default -1. end of hold length
        if endhole>-1 a gap will be added if TO.last< endhole +/- tol
        if endhole>-1 and TO.last> endhole +/- tol a warning will be raised
    
    Returns:
    np_nf,np_nt,np_nID: numpy 1D arrays with new from, to and interval ID
    np_gap,np_overlap: numpy 1D arrays with FROM position where gaps and 
                       overlaps where detected.
    
    """
    
    cdef: 
        int i
        int nint=-1
        int ngap=-1
        int noverlap=-1
        int ndat= len(in_f)
    
    
    #make a deep copy of the from to intervals and create a memory view
    np_f= np.zeros(ndat, dtype=float)
    np_t= np.zeros(ndat, dtype=float)
    cdef:
        double[:] f = np_f
        double[:] t = np_t
    for i in range(ndat):
        f[i]=in_f[i]
        t[i]=in_t[i]
    
    #make a long array (reserve memory space)
    np_nf= np.zeros(ndat*2+4, dtype=float)
    np_nt= np.zeros(ndat*2+4, dtype=float)
    np_nID= np.zeros(ndat*2+4, dtype=int)
    np_gap= np.zeros(ndat*2+4, dtype=int)
    np_overlap= np.zeros(ndat*2+4, dtype=int)
    cdef:
        double[:] nt = np_nt
        long[:] nID = np_nID
        long[:] gap = np_gap
        double[:] nf = np_nf
        long[:] overlap = np_overlap
    
    # gap first interval 
    if f[0]>tol:
        nint+=1
        nf[nint]=0.
        nt[nint]=f[0]
        nID[nint]=-999
        ngap+=1
        gap[ngap]=id[0]
               
    
    for i in range(ndat-1):
        
        # there is no gap?
        if -tol<=f[i+1]-t[i]<=tol: 
            # existing sample
            nint+=1
            nf[nint]=f[i]
            nt[nint]=f[i+1]
            nID[nint]=id[i]
            continue

        # there is a gap?
        if f[i+1]-t[i]>=tol: 
            # existing sample
            nint+=1
            nf[nint]=f[i]
            nt[nint]=t[i]
            nID[nint]=id[i]
            #gap
            nint+=1
            nf[nint]=t[i]
            nt[nint]=f[i+1]
            nID[nint]=-999
            ngap+=1
            gap[ngap]=id[i]
            continue 
            
        # there is overlap?
        if f[i+1]-t[i]<=-tol: 
            # the overlap is smaller that the actual sample? 
            if f[i+1]>f[i]: 
                # existing sample
                nint+=1
                nt[nint]=max(f[i+1],f[i]) # ising max here to avoid negative interval from>to
                nf[nint]=f[i]
                nID[nint]=id[i]
                noverlap+=1
                overlap[noverlap]=id[i]
                continue 
            # large overlap?
            else: 
                #whe discard next interval by making it 0 length interval
                # this will happend only in unsorted array... 
                noverlap+=1
                overlap[noverlap]=id[i]
                # update to keep consistency in next loopondh
                nint+=1
                nt[nint]=t[i] # ising max here to avoid negative interval from>to
                nf[nint]=f[i]
                nID[nint]=id[i]
                f[i+1]=t[i+1]
               
    # add last interval
    # one sample drillhole? 
    '''
    if ndat==1:
        nint+=1
        nt[nint]=t[0] # ising max here to avoid negative interval from>to
        nf[nint]=f[0]
        nID[nint]=id[0]
    '''
        
    # there are not problem (like a gap or an overlap)
    if (-tol<=f[ndat-1]-t[ndat-2]<=tol) and ndat>1:
        nint+=1
        nt[nint]=t[ndat-1] # ising max here to avoid negative interval from>to
        nf[nint]=f[ndat-1]
        nID[nint]=id[ndat-1]
    else:
        # just add the sample (the problem was fixed in the previous sample)
        nint+=1
        nt[nint]=t[ndat-1] # ising max here to avoid negative interval from>to
        nf[nint]=f[ndat-1]
        nID[nint]=id[ndat-1]
    
    # add end of hole
    if endhole>-1:
        # there is an end of hole gap?
        if (tol>endhole-t[ndat-2]) and ndat>1:
            nint+=1
            nt[nint]=endhole
            nf[nint]=t[ndat-1]
            nID[nint]=-999
            ngap+=1
            gap[ngap]=-888  # this is a gap at end of hole
        # there is an end of hole overlap?
        if (tol>endhole-t[ndat-2]) and ndat>1:
            nint+=1
            nt[nint]=endhole
            nf[nint]=t[ndat-1]
            nID[nint]=-999
            noverlap+=1
            overlap[noverlap]=-888 # this is an overlap at end of hole
        
        # there is no gap or overlap, good... then fix small differences
        if (tol<endhole-t[ndat-2]) and (endhole-t[ndat-2]>-tol):
            nt[nint]=endhole

    # make first interval start at zero, if it is != to zero but close to tolerance
    if 0<nf[0]<=tol:
        nf[0]=0
        


        
    
    return np_nf[:nint+1],np_nt[:nint+1],np_nID[:nint+1],np_gap[:ngap+1],np_overlap[:noverlap+1]



#-------------------------------------------------------------------
#  Drillhole class
#-------------------------------------------------------------------
cdef class Drillhole:
    """
    Drillhole(collar, survey)
    
    Drillhole working database object with functions to desurvey and
    validate drillholes.
    
    Parameters
    ----------
    collar : Pandas DataFrame 
        Collar table containing compulsory fields: BHID with any dtype
        and XCOLLAR,YCOLLAR, ZCOLLAR with dtypes float64.  
        Collar table containing compulsory field: LENGTH with dtypes 
        float64. This is the length of drillhole and can be used in 
        some functions, for example, to fill gaps. 
    survey : Pandas DataFrame 
        Survey table containing compulsory fields: BHID with any dtype
        and AT,AZ, DIP with dtypes float64. 
    
    Attributes
    ----------
    collar : Pandas DataFrame 
    survey : Pandas DataFrame 
    tables : [Pandas DataFrame]
        list with Pandas DataFrame with interval tables (ex. assay)
        containing compulsory fields: BHID with any dtype
        and FROM, TO with dtypes float64.
    table_mames : [str]
        list of table names
    
    Notes
    -----
    A new copy of the input data will be created in memory. To work with 
    shared memory in an external DataFrame you may copy back the table:
    
    ex. >> shared_collar = mydrillholeDB.collar 
    
    1) To add interval tables use addtable
    2) Only the existence of compulsory fields are validated in 
       object initialization and when adding interval tables
    3) Survey tables may have at least 2 interval and interval AT=0 is 
       required for desurvey
    
    """ 
    cdef readonly object collar
    cdef readonly object survey
    cdef readonly object table
    
    property table_mames:
        def __get__(self):
            return self.table.keys()
    
    def __cinit__(self, collar, survey):
        """
        add doc string here
        """
        #check the input is correct
        assert isinstance(collar, pd.DataFrame), "collar is not a pandas DataFrame" 
        assert isinstance(survey, pd.DataFrame), "survey is not a pandas DataFrame" 
        
        #check we have the rigth naming in collar columns 
        assert 'BHID' in collar.columns, "collar don't have BHID column"
        assert 'XCOLLAR' in collar.columns, "collar don't have XCOLLAR column"
        assert 'YCOLLAR' in collar.columns, "collar don't have YCOLLAR column"
        assert 'ZCOLLAR' in collar.columns, "collar don't have ZCOLLAR column"
        
        
        if 'LENGTH' not in collar.columns: 
            warnings.warn('! Collar table without LENGTH field' ) 
        
        #check we have the rigth naming in survey columns 
        assert 'BHID' in survey.columns, "survey don't have BHID column"
        assert 'AT' in survey.columns, "survey don't have AT column"
        assert 'AZ' in survey.columns, "survey don't have AZ column"
        assert 'DIP' in survey.columns, "survey don't have DIP column"
        
        
        self.collar = collar.copy(deep=True) # we remove external reference
        self.survey = survey.copy(deep=True) # we remove external reference
        self.table = {}
        
        # set numeric columns as float (we don't want integers)
        #this rize an error if there are non numeric values
        self.collar['XCOLLAR'] = self.collar['XCOLLAR'].astype(float)
        self.collar['YCOLLAR'] = self.collar['YCOLLAR'].astype(float)
        self.collar['ZCOLLAR'] = self.collar['ZCOLLAR'].astype(float)
        if 'LENGTH' in self.collar.columns: 
            self.collar['LENGTH'] = self.collar['LENGTH'].astype(float)
        self.survey['AT'] = self.survey['AT'].astype(float)
        self.survey['AZ'] = self.survey['AZ'].astype(float)
        self.survey['DIP'] = self.survey['DIP'].astype(float)
        
        
        
        
        # set BHID as uppercase to avoid issues
        self.collar['BHID']= self.collar['BHID'].str.upper()
        self.survey['BHID']= self.survey['BHID'].str.upper()
        
        # sort the data 
        self.collar.sort_values(by=['BHID'], inplace=True)
        self.survey.sort_values(by=['BHID','AT'], inplace=True)
           
    
    cpdef addtable(self,object table,str table_name,bint overwrite =False):
        """
        addtable(table, table_name,overwrite = False)
        
        Add an interval table and assign a name.
        
        Parameters
        ----------
        table : Pandas DataFrame 
                containing compulsory fields BHID with any dtype and 
                FROM, TO with dtypes float64.
        table_name : str 
                with an the name of the table
        overwrite : boolean, optional, default False
                if the name of the table exists overwrite = True will 
                overwrite the existing table with the new input table
                
        Examples
        --------
        
        >>> mydrillhole.addtable(assay, 'assay')
        
                
        """        
        #check the input is correct
        assert isinstance(table, pd.DataFrame), "table is not a pandas DataFrame"
        #assert isinstance(table_name, str), "table_name is not a string"
        
        #check we have the rigth naming in columns 
        assert 'BHID' in table.columns, "%s don't have BHID" %table_name
        assert 'FROM' in table.columns, "%s don't have FROM" %table_name
        assert 'TO' in table.columns, "%s don't have TO" %table_name
        
        
        if table_name not in self.table:
            self.table[table_name]=table.copy(deep=True) # we remove external reference
        else:
            if overwrite == True: 
                self.table[table_name]=table.copy(deep=True) # we remove external reference
            else:
                raise NameError('Table %s already exist, use overwrite = True to overwrite' % table_name)


        # set numeric columns as float (we don't want integers)
        #this rize an error if there are non numeric values
        self.table[table_name]['FROM'] = self.table[table_name]['FROM'].astype(float)
        self.table[table_name]['TO'] = self.table[table_name]['TO'].astype(float)

            
        # set BHID as uppercase to avoid issues
        self.table[table_name]['BHID']= self.table[table_name]['BHID'].str.upper()
        
        # sort the data 
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)
        
        # reset index, if necessary uncomment next line
        # self.table[table_name].reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')

    cpdef del_table(self,str table_name):
        """
        del_table(table_name)
        
        Delete an interval table from the drillhole database object.
        
        Parameters
        ----------
        table_name : str 
                with an the name of the table
                
        Examples
        --------
        
        >>> mydrillhole.addtable('assay')
        
                
        """        
        #check the input is correct
        if table_name in self.table:
            del self.table[table_name]
        else:
            raise NameError('Table {} not found in the drillhole database object' % table_name)
            


    cpdef validate(self):
        """
        validate()
        
        Runs a set of basic validations on survey and collar 
        consisting in: 
        - check existence of Null values
        - check survey without values AT=0
        - check that coordinates and direction angles dtypes are float64
        - check survey without collar and collar without survey 
        
        Notes
        -----
        a) You may run this validation before doing desurvey
        b) Only few minimum validation are implemented for now
        c) No value is returned, it raises an error 
        
        TODO
        ----
        Implement check relation between table, large variations in
        survey angles, missing BHID.
        
        Collect all errors and return a list of errors. 
        Collect all warnings and return a list of warnings. 
        
        See Also
        --------
        validate_table
        
        Examples
        --------
        
        >>> mydrillhole.validate()
        
                
        """  
        
        # Warning:  hasnans() is for pandas 0.16, hasnans for pandas 0.17
        
        #check collar
        
        #check repeated collar BHID
        if len(self.collar.loc[self.collar.duplicated(['BHID'])])>0:
            raise NameError('There are duplicated BHID at collar')

        #check repeated survey BHID,AT
        if len(self.survey.loc[self.survey.duplicated(['BHID','AT'])])>0:
            raise NameError('There are duplicated BHID,AT at survey')
        
        # null values in collar
        if self.collar['BHID'].hasnans:
            raise NameError('Non defined BHID in collar table')
        if self.collar['XCOLLAR'].hasnans:
            raise NameError('Non defined XCOLLAR in collar table')
        if self.collar['YCOLLAR'].hasnans:
            raise NameError('Non defined YCOLLAR in collar table')
        if self.collar['ZCOLLAR'].hasnans:
            raise NameError('Non defined ZCOLLAR in collar table')
        if self.collar['XCOLLAR'].dtypes!='float64':
            raise NameError('XCOLLAR in collar table != float64')            
        if self.collar['YCOLLAR'].dtypes!='float64':
            raise NameError('YCOLLAR in collar table != float64')
        if self.collar['ZCOLLAR'].dtypes!='float64':
            raise NameError('ZCOLLAR in collar table != float64')
        
        if 'LENGTH' in self.collar.columns:         
            if self.collar['LENGTH'].dtypes!='float64':
                raise NameError('LENGTH in collar table != float64') 
            if self.collar['LENGTH'].hasnans:
                raise NameError('Non defined LENGTH values in collar table')
            

        #check SURVEY
        # null values in survey   
        
        if self.survey['BHID'].hasnans:
            raise NameError('Non defined BHID in survey table')
        if self.survey['AT'].hasnans:
            raise NameError('Non defined AT in survey table')
        if self.survey['AZ'].hasnans:
            raise NameError('Non defined AZ in survey table')
        if self.survey['DIP'].hasnans:
            raise NameError('Non defined DIP in survey table')
        
        #check dtypes 
        if self.survey['AT'].dtypes!='float64':
            raise NameError('AT in survey table != float64')            
        if self.survey['DIP'].dtypes!='float64':
            raise NameError('DIP in survey table != float64')
        if self.survey['AZ'].dtypes!='float64':
            raise NameError('AZ in survey table != float64')
        
        #Check using same data type for BHID
        if self.survey['BHID'].dtypes!=self.collar['BHID'].dtypes:
            raise NameError("self.survey['BHID'].dtypes!=self.collar['BHID'].dtypes")
        
        #check survey without values: at=0
        self.survey.sort_values(by=['BHID','AT'], inplace=True)
        error = self.__checkAt0(self.survey['BHID'].values, self.survey['AT'].values)
        if error>-1:
            raise NameError('Firts inteval AT!=0 at survey table, positiom %d' %error) 
        
        #check survey with only one survey occurrence per drillholes: this produces error in desurvey 
        for i in self.survey['BHID'].unique():
            ns=self.survey['BHID'][self.survey['BHID'] == i].count()
            if ns<2:
                warnings.warn('! survey with one value at BHID: {}. This will produce error at desurvey'.format(i) )
        
        
        #check survey without collar
        cID = self.collar['BHID'].values
        for i in self.survey['BHID'].unique(): 
            if i not in cID:
                warnings.warn('! survey without collar at BHID: {}'.format(i) )    
        
        #check collar without survey
        sID = self.survey['BHID'].unique()
        for i in self.collar['BHID'].values: 
            if i not in sID:
                warnings.warn('! collar without survey at BHID: {}'.format(i) )    
    
            
        return None
        # TODO: check table relationship
        # TODO: check survey.AT.last > endofhole
        
    cdef __checkAt0(self,np.ndarray BHID, np.ndarray[double, ndim=1] AT):
        # this is a hide function
        # the input data is assumed sorted
        cdef int n= AT.shape[0]
        cdef int i, start
    
        
        # first interval is zero
        if  AT[0]>0.00001:
            return 0
            
        # the first dhole intervals (where BHID[i-1]!=BHID[i]) are zero?
        for i in range(1,n):
            if BHID[i-1]!=BHID[i]:
                if  AT[i]>0.00001: 
                   return i
        
        return -1            

    # TODO: Optimize this one, it is too slow
    cpdef fix_survey_one_interval_err(self, double dummy_at):
        """
        fix_survey_one_interval_err(dummy_at)
        
        This function add a dummy survey in records with one survey only
        
        The desurvey algorithm may produce unexpected results if there is only 
        one survey interval per drillhole. This function avoid this issue
        by duplication the survey parameter at a dummy_at position. 
        
        TODO: use collar length instead dummy_at, if LENGTH is defined
        
        """
        
        
        #check survey with only one survey occurrence per drillholes: 
        #this produces error in desurvey 
        for i in self.survey['BHID'].unique():
            ns=self.survey['BHID'][self.survey['BHID'] == i].count()
            bhid=[]
            at = []
            az =[]
            dip =[]
            if ns==1:
                bhid.append(self.survey['BHID'][self.survey['BHID'] == i].iloc[0])
                az.append(self.survey['AZ'][self.survey['BHID'] == i].iloc[0])
                dip.append(self.survey['DIP'][self.survey['BHID'] == i].iloc[0])
                at.append(dummy_at)
            
            tmp=pd.DataFrame({'BHID':bhid, 'AZ': az, 'DIP': dip, 'AT': at })
            
            # append this to survey 
            self.survey=self.survey.append(tmp)
            
            #sort survey
            self.survey.sort_values(by=['BHID','AT'], inplace=True)
                
        
    cpdef validate_table(self, table_name): 
        """
        validate_table()
        
        Runs a set of basic validations on interval tables 
        consisting in: 
        - check null values in table BHID
        - check null values in From/To
        - check that FROM and TO dtypes are float64
        
        Notes
        -----
        a) You may run this validation before doing desurvey
        b) Ony few minimum validation are implemented for now
        c) No value is returned, it raises an error 
        
        TODO
        ----
        Implement check relation between tables, gaps, overlays, 
        missing BHID, missing Survey, Interval longer that max_depth.
        
        Collect all errors and return a list of errors. 
        Collect all warnings and return a list of warnings. 
        
        See Also
        --------
        validate
                
        Examples
        --------
        
        >>> mydrillhole.validate_table('assay')
        
                
        """  
           
        #check the input is correct
        assert isinstance(table_name, str), 'table_name is not a string'
        assert table_name in self.table, '%s not exist in this drillhole database' % table_name
        
        #check table

        #check repeated BHID,FROM
        if len(self.table[table_name][self.table[table_name].duplicated(['BHID','FROM'])])>0:
            raise NameError('There are duplicated BHID,FROM at table')


        # null values in table bhid
        if self.table[table_name]['BHID'].hasnans:
            raise NameError('Non defined BHID in %s' % table_name)
        # null values in From/To
        if self.table[table_name]['FROM'].hasnans:
            raise NameError('Non defined FROM in %s' % table_name)
        if self.table[table_name]['TO'].hasnans:
            raise NameError('Non defined TO in %s' % table_name)
        if self.table[table_name]['FROM'].dtypes!='float64':
            raise NameError('FROM in table %s != float64' % table_name)
        if self.table[table_name]['TO'].dtypes!='float64':
            raise NameError('TO in table %s != float64' % table_name)

        #Check using same data type for BHID
        if self.table[table_name]['BHID'].dtypes!=self.collar['BHID'].dtypes:
            raise NameError("self.table[%s]['BHID'].dtypes!=self.collar['BHID'].dtypes" % table_name)

        #check table without collar
        cID = self.collar['BHID'].values
        for i in self.table[table_name]['BHID'].unique(): 
            if i not in cID:
                warnings.warn('Warning on table {}: BHID: {} not in collar'.format(table_name, i) )    
        
        #check collar without table
        tID = self.table[table_name]['BHID'].unique()
        for i in self.collar['BHID'].values: 
            if i not in tID:
                warnings.warn('Info : collar BHID {} not at table {}'.format(i,table_name) )    

        
        
        # TODO: check overlaps and table relationship
        # TODO: check TO.last > endofhole

    
    cpdef txt2intID(self, str table_name):
        """
        txt2intID(table_name)
        
        Creates an alternative BHID of type integer on the 
        table[table_name] and in collar. The new field BHIDint is just
        and ordered list of integers. 
        
        BHIDint may be required in some functions compiles in fortran 
        for pyslib.gamv.
        
        Parameters
        ----------
        table_name : str
        
        
        Notes
        -----
        The Collar and the table may be sorted. 
        
        Examples
        --------
        
        >>> mydrillhole.collar.sort(['BHID'], inplace=True)
        >>> mydrillhole.table['assay'].sort(['BHID', 'FROM'], inplace=True)
        >>> mydrillhole.txt2intID('assay')
        
        """

        # the data is assumed sorted
        
        assert table_name in self.table, 'The table %s do not exist in the database' % table_name
        
        cdef int sloc
        cdef int i
        cdef int j
        cdef int nc = self.collar.shape[0]
        cdef int nt = self.table[table_name].shape[0]
        cdef np.ndarray[long, ndim=1] cBHID = np.zeros([nc], dtype=long)
        cdef np.ndarray[long, ndim=1] tBHID = np.zeros([nt], dtype=long)
        cdef np.ndarray[object, ndim=1] cTextID = np.empty([nc], dtype=object)
        cdef np.ndarray[object, ndim=1] tTextID = np.empty([nt], dtype=object)
        
        cTextID[:] = self.collar['BHID']
        tTextID[:] = self.table[table_name]['BHID']
        

        sloc = 0
        for i in range(nc):

            cBHID[i] = i+1

            # get first position
            for j in range(sloc, nt): 
                if cTextID[i] != tTextID[j]:
                    continue
                else:
                    sloc = j
                    break
            # operate in the first position 
            for j in range(sloc, nt): 
                if cTextID[i] == tTextID[j]:
                    tBHID[j]=cBHID[i] 
                else:
                    sloc = j
                    break  
        
        self.collar['BHIDint']= cBHID
        self.table[table_name]['BHIDint'] = tBHID



    cpdef desurvey(self, str table_name, bint endpoints=False, 
                   bint warns=True):
        """
        desurvey(table_name, endpoints=False)
        
        Create coordinates and direction angles at table midpoint 
        intervals. If endpoints=True it also create coordinate fields
        at end point intervals
        
        This function calls the function desurv1dh at each interval 
        point of the table. The result are added as new fields. If the 
        coordinate field exist they will be overwrited
        
        
        Parameters
        ----------
        table_name : str
        endpoints : boolean 
        
        See Also
        --------
        dsmincurb, desurv1dh
        
        
        Notes
        -----
        The Collar, Survey and the table may be sorted. 
        If you call the function with endpoints=False end points already
        desurveyed may not be overwrited. 
        
        Examples
        --------
        
        >>> mydrillhole.collar.sort(['BHID'], inplace=True)
        >>> mydrillhole.survey.sort(['BHID', 'AT'], inplace=True)
        >>> mydrillhole.table['assay'].sort(['BHID', 'FROM'], inplace=True)
        >>> mydrillhole.desurvey('assay', endpoints=True)
        
        """
        
        
        #check the input is correct
        assert table_name in self.table, "table %s not exist" % table_name
        
        
        cdef idc = self.collar['BHID'].values
        cdef np.ndarray[double, ndim=1] xc = self.collar['XCOLLAR'].values
        cdef np.ndarray[double, ndim=1] yc = self.collar['YCOLLAR'].values
        cdef np.ndarray[double, ndim=1] zc = self.collar['ZCOLLAR'].values
        cdef ids = self.survey['BHID'].values
        cdef np.ndarray[double, ndim=1] ats = self.survey['AT'].values
        cdef np.ndarray[double, ndim=1] azs = self.survey['AZ'].values
        cdef np.ndarray[double, ndim=1] dips = self.survey['DIP'].values
        cdef idt =self.table[table_name]['BHID'].values
        cdef np.ndarray[double, ndim=1] fromt = self.table[table_name]['FROM'].values
        cdef np.ndarray[double, ndim=1] tot = self.table[table_name]['TO'].values 
                         
        # internal
        cdef int nc= idc.shape[0]
        cdef int ns= ids.shape[0]
        cdef int nt= idt.shape[0]
        cdef int jc
        cdef int js
        cdef int jt
        cdef int indbs
        cdef int indbt
        cdef int indes
        cdef int indet
        cdef int inds,
        cdef int indt
        cdef float mid
        
        # otput
        cdef np.ndarray[double, ndim=1] azmt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] dipmt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] xmt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] ymt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] zmt = np.empty([nt], dtype=float)
                
        
        #if endpoints==true:
        cdef float tmpaz, tmpdip
        cdef np.ndarray[double, ndim=1] xbt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] ybt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] zbt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] xet = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] yet = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] zet = np.empty([nt], dtype=float)
        
        # initialize output to nans
        azmt[:] = np.nan
        dipmt[:] = np.nan
        xmt[:] = np.nan
        ymt[:] = np.nan
        zmt[:] = np.nan
        xbt[:] = np.nan
        ybt[:] = np.nan
        zbt[:] = np.nan
        xet[:] = np.nan
        yet[:] = np.nan
        zet[:] = np.nan
              
        indbt = 0
        indet = 0 
        inds = 0
        indt = 0   
        for jc in range(nc): # for each collar 
            indbs = -1
            indes = -1
            # get first index of the collar jc in survey 
            for js in range(inds, ns):
                # find index begin and end for the actual collar
                if idc[jc]==ids[js]:
                    inds = js
                    indbs = js
                    break
            
            # get last index of the collar jc in survey 
            for js in range(inds, ns):
                # find index begin and end for the actual collar
                if idc[jc]!=ids[js]:
                    break
                else: 
                    inds = js
                    indes = js
            
            # this is not working but check included in validation 
            if indbs==-1 or indes==-1:
                # do not desurvey this drillhole
                warnings.warn('! collar {} without survey, table not desurveyed'.format(idc[jc]))
                continue

            # this is not working but check included in validation 
            if indbs==indes:
                # do not desurvey this drillhole
                warnings.warn('! collar {} without survey at end collar, table not desurveyed'.format(idc[jc]))
                continue

            
            # with the index indbs and indes we desurvey each collar
            for jt in range(indt, nt):
                # the table id is similar to collar? Then desurvey
                
                if idc[jc]==idt[jt]:
                    #desurvey this point
                    indt = jt # do not loop again before this index
                    
                    # desurvey mid interv                    
                    mid = fromt[jt] + (tot[jt]-fromt[jt])/2.
                    
                    azmt[jt],dipmt[jt],xmt[jt],ymt[jt],zmt[jt] = \
                    desurv1dh(indbs,indes,ats,azs,dips,xc[jc],yc[jc],zc[jc],mid,warns)
                    
                    if endpoints==True:
                        tmpaz,tmpdip,xbt[jt],ybt[jt],zbt[jt] = \
                        desurv1dh(indbs,indes,ats,azs,dips,xc[jc],yc[jc],zc[jc],fromt[jt],warns)
                        
                        tmpaz,tmpdip,xet[jt],yet[jt],zet[jt] = \
                        desurv1dh(indbs,indes,ats,azs,dips,xc[jc],yc[jc],zc[jc],tot[jt],warns)  
                        
        
        self.table[table_name]['azm'] = azmt
        self.table[table_name]['dipm']= dipmt
        self.table[table_name]['xm']= xmt
        self.table[table_name]['ym']= ymt
        self.table[table_name]['zm']= zmt
        if endpoints==True:
            self.table[table_name]['xb']= xbt
            self.table[table_name]['yb']= ybt
            self.table[table_name]['zb']= zbt
            self.table[table_name]['xe']= xet
            self.table[table_name]['ye']= yet
            self.table[table_name]['ze']= zet


        # this is to produce a warning if some intervals where not desurvey
        try: 
            assert  np.isfinite(self.table[table_name]['xm'].values).all(), "no finite coordinates at xb, please filter out non finite coordinates and try again"
        except:
            warnings.warn('Some intervals where non-desurveyed and NAN coordinates where created, check errors in collar or survey')



    #-------------------------------------------------------------------
    #       Table gap, compositing and merging
    #-------------------------------------------------------------------
    cpdef add_gaps(self,  str table_name, 
                          str new_table_name, 
                          bint overwrite =False,
                          double tol=0.01, 
                          bint endhole=False,
                          bint clean=True):
        
        """
        mydrillhole.add_gaps(table_name, new_table_name,
                             overwrite =False, tol=0.01,
                             endhole=False, clean=True)
        
        Fill gaps in one drillhole with new FROM-TO intervals, including 
        the gaps at collar and at the end of drillholes. 
        
        _id0 will be added to the new table and the existing table, 
        it links input table rows with the output table rows, gaps 
        have _id0= -999
        
               
        Input:
        table_name: name of the table to add gaps
            it must be an existing table at drillhole.table.keys()
        new_table_name: name of the new table with gaps added
            it can be a new name, no in drillhole.table.keys()
            or a table name on drillhole.table.keys() if overwrite =True
        overwrite: default True. 
            If the table exist and overwrite = True the existing table 
            will be overwrite. 
        tol: default 0.01. 
            gaps and overlays within tol will be ignore but adjusted
            in order to impose FROM[i-1]==TO[i]
        endhole: default False.
            if true a gap at the end of the drillhole will be added. 
            The end of the drillhole is calculated from 
            drillhole.collar['LENGTH'], and error will be raised if:
            there is no 'LENGTH' field at collar or if there are 
            undefined values at LENGTH. Warnings will be raised is the 
            TO  > LENGTH.
        clean: default True.
            Delete temporary columns created with suffix __tmp__. 
        
        Returns:
        
        gap,overlap: _id0 of the raw before the gap or the overlap was 
                     detected. Gaps and overlaps ate the end of the 
                     drillhole (if endhole==True) have value -888.
        
        """
        
        cdef:
            double l_endhole
            
            
        
        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)
        
        if endhole==True: 
            assert 'LENGTH' in self.collar.columns, 'There is no LENGTH field at collar, use endhole=False'
        
        # sort table
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)
        
        # add ID to the table
        self.table[table_name].loc[:,'_id0']= np.arange(self.table[table_name].shape[0])[:]
        
        # create a group to easily iterate
        group=self.table[table_name].groupby('BHID')
        
        #add gaps
        BHID=group.groups.keys()
        nnf=[]
        nnt=[]
        nnID=[]
        nnBHID=[]
        nngap= []
        nnoverlap = []
        for i in BHID:
            if endhole:
                l_endhole=group.get_group(i)['LENGTH'].values
            else:
                l_endhole=-1
            
            nf,nt,nID,gap,overlap=fillgap1Dhole(in_f = group.get_group(i)['FROM'].values, 
                                          in_t = group.get_group(i)['TO'].values, 
                                          id = group.get_group(i)['_id0'].values, 
                                          tol=tol,
                                          endhole=l_endhole)
            
            
            nBHID = np.empty([len(nf)], dtype=object, order='C')
            nBHID[:]=i
            nnf+=nf.tolist()
            nnt+=nt.tolist()
            nnID+=nID.tolist()
            nnBHID+=nBHID.tolist()
            nngap+=gap.tolist()
            nnoverlap+=overlap.tolist()
    
    
        #create new table with gaps (only with fields )
        newtable=pd.DataFrame({'BHID':nnBHID, 'FROM':nnf,'TO':nnt,'_id0':nnID})
        
        newtable=newtable.join(self.table[table_name], on='_id0', rsuffix='__tmp__')
        
        #clean if necessary
        if clean:
            newtable.drop(
               ['BHID__tmp__', 'FROM__tmp__','TO__tmp__','_id0__tmp__'],
               axis=1,inplace=True, errors='ignore')
        
        #add table to the class
        self.addtable(newtable,new_table_name,overwrite)
                
        return nngap,nnoverlap


    cpdef merge(self,str table_A,
                     str table_B,
                     str new_table_name,
                     bint overwrite =False, 
                     double tol=0.01,
                     bint clean=True):
        
        """
        mydrillhole.merge(table_A, table_B, new_table_name,
                          overwrite =False, tol=0.01, clean=True)
        
        Combine two tables in one by intersecting intervals. 
        
        This function requires drillholes without gaps and overlaps. 
        run add_gaps in table_A and table_B before using this function
        
               
        Input:
        table_A: name of the first table
            it must be an existing table at drillhole.table.keys()
        table_B: name of the second table
            it must be an existing table at drillhole.table.keys()
        new_table_name: name of the new table
            it may not exists at drillhole.table.keys()
        overwrite: default True. 
            If new_table_name exists and overwrite == True the existing 
            table will be overwrite. 
        tol: default 0.01. 
            segments, gaps and overlays within tol will be ignore but 
            adjusted
        clean: default True.
            Delete temporary columns created with suffix __tmp__. 
        
        Returns:
        
        gap,overlap: _id0 of the raw before the gap or the overlap was 
                     detected. Gaps and overlaps ate the end of the 
                     drillhole (if endhole==True) have value -888.
        
        """
        
        
        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)

        
        # sort tables
        self.table[table_A].sort_values(by=['BHID', 'FROM'], inplace=True)
        self.table[table_B].sort_values(by=['BHID', 'FROM'], inplace=True)
        
        # add ID to tables
        self.table[table_A].loc[:,'_id0']= np.arange(self.table[table_A].shape[0])[:]
        self.table[table_B].loc[:,'_id1']= np.arange(self.table[table_B].shape[0])[:]
        
        # create a groups to easily iterate
        groupA=self.table[table_A].groupby('BHID')
        groupB=self.table[table_B].groupby('BHID')
        
        
        # prepare fixed long array to send data    
        #    input     
        np_la=np.empty(self.table[table_A].shape[0], dtype = float)
        np_lb=np.empty(self.table[table_B].shape[0], dtype = float)
        np_ida=np.empty(self.table[table_A].shape[0], dtype = int)
        np_idb=np.empty(self.table[table_B].shape[0], dtype = int)
        
        ll = self.table[table_A].shape[0] +  self.table[table_B].shape[0] +10
        
        nBHID = np.empty(ll, dtype=object, order='C')
        
        cdef:
            double[::1] la = np_la
            double[::1] lb = np_lb
            long[::1] ida  = np_ida
            long[::1] idb  = np_idb
            double[::1] lab
            long[::1] newida
            long[::1] newidb 
            
            int k
            int nk

            int j
            int nj
            
            int na
            int nb
            
            int n
            
            double endf
        
        #merge
        BHID=self.collar.BHID.values
        nnf=[]
        nnt=[]
        nnIDA=[]
        nnIDB=[]
        nnBHID=[]
        
        keysA= groupA.groups.keys()
        keysB= groupB.groups.keys()
        
        for i in BHID:
            
            
            # if we really have to merge
            if (i in keysA) and (i in keysB): 
            
                # prepare input data
                # table A drillhole i
                nk=groupA.get_group(i).shape[0]
                for k in range(nk):
                    la[k]=groupA.get_group(i)['FROM'].values[k]
                    ida[k]=groupA.get_group(i)['_id0'].values[k]
                
                la[nk]=groupA.get_group(i)['TO'].values[nk-1]
                ida[nk]=groupA.get_group(i)['_id0'].values[nk-1]
                
                # table B drillhole i
                nj=groupB.get_group(i).shape[0]
                for j in range(nj):
                    lb[j]=groupB.get_group(i)['FROM'].values[j]
                    idb[j]=groupB.get_group(i)['_id1'].values[j]
                
                lb[nj]=groupB.get_group(i)['TO'].values[nj-1]
                idb[nj]=groupB.get_group(i)['_id1'].values[nj-1]
                
                # make sure the two drill holes have the same length 
                # by adding a gap at the end of the shortest drillhole
                if lb[nj] > la[nk]:
                    nk+=1
                    la[nk] = lb[nj]
                    ida[nk] = -999
                    endf = lb[nj]
                
                elif la[nk] > lb[nj]:
                    nj+=1
                    lb[nj] = la[nk]
                    idb[nj] = -999
                    endf = la[nk]
                
                # merge drillhole i
                n, np_lab, np_newida, np_newidb = merge_one_dhole(la[:nk+1],lb[:nj+1], ida[:nk+1], idb[:nj+1], tol=0.01)

                # dhid
                nBHID[:n]=i
                nnBHID+=nBHID[:n].tolist()
                
                # from
                nnf+=np_lab[:-1].tolist()
                # to
                nnt+=np_lab[1:].tolist() 
                
                # IDs
                nnIDA+=np_newida[:-1].tolist()
                nnIDB+=np_newidb[:-1].tolist()
            
                continue
            
            # it is only on table A?
            if (i in keysA): 
                
                n= groupA.get_group(i).shape[0]
                
                # in this case we add table A and ignore B
                # dhid
                
                nBHID[:n]=i
                nnBHID+=nBHID[:n].tolist()
                
                # from
                nnf+=groupA.get_group(i)['FROM'].values.tolist()
                # to
                nnt+=groupA.get_group(i)['TO'].values.tolist() 
                
                # IDs
                tmp=-999*np.ones(n, dtype='int')
                nnIDA+=groupA.get_group(i)['_id0'].values.tolist()
                nnIDB+= tmp.tolist() 
                continue
            
            # it is only on table B?
            if (i in keysB): 
                
                n= groupB.get_group(i).shape[0]
                
                # in this case we add table B and ignore A
                # dhid
                
                nBHID[:n]=i
                nnBHID+=nBHID[:n].tolist()
                
                # from
                nnf+=groupB.get_group(i)['FROM'].values.tolist()
                # to
                nnt+=groupB.get_group(i)['TO'].values.tolist() 
                
                # IDs
                tmp=-999*np.ones(n, dtype='int')
                nnIDA+= tmp.tolist() 
                nnIDB+= groupB.get_group(i)['_id1'].values.tolist() 
                continue
                

        
    
        #create new table with intervals and ID
        newtable=pd.DataFrame({'BHID':nnBHID, 'FROM':nnf,'TO':nnt,'_id0':nnIDA,'_id1':nnIDB})
        
        # merge with existing data
        newtable=newtable.join(self.table[table_A], on='_id0', rsuffix='__tmp__')
        newtable=newtable.join(self.table[table_B], on='_id1', rsuffix='__tmp__')

        
        #clean if necessary
        if clean:
            newtable.drop(
               ['BHID__tmp__', 'FROM__tmp__','TO__tmp__','_id0__tmp__','_id1__tmp__'],
               axis=1,inplace=True, errors='ignore')
        
        #add table to the class
        self.addtable(newtable,new_table_name,overwrite)


    cpdef fix_zero_interval(self, str table_name, 
                          str new_table_name, 
                          bint overwrite = False,
                          double tol = 0.01,
                          bint clean = True,
                          bint endhole=False,
                          bint addgaps=True):
        
        
        """
        gap_assay,overlap_assay = mydrillhole.fix_zero_interval(
                          table_name, new_table_name,
                          overwrite =False, tol=0.01, clean = True,
                          endhole=False, addgaps=True)
        
        Remove zero length intervals 
        
        The function first removes zero intervals and then  
        add gaps if the gap length is longer than the tolerance. 
        
        Note: 
        This function remove zero intervals and the calls the add_gaps
        function. 
        
        Input:
        table_name: name of the table
        
        new_table_name: name of the new table
            it may not exists at drillhole.table.keys()
        overwrite: default True. 
            If new_table_name exists and overwrite == True the existing 
            table will be overwrite. 
        tol: default 0.01. 
            segments with length<= 0.01 will be considered as zero
            length and removed
        endhole: default False.
            if true a gap at the end of the drillhole will be added. 
            The end of the drillhole is calculated from 
            drillhole.collar['LENGTH'], and error will be raised if:
            there is no 'LENGTH' field at collar or if there are 
            undefined values at LENGTH. Warnings will be raised is the 
            TO  > LENGTH.
        clean: default True.
            Delete temporary columns created with suffix __tmp__.                         
        """

        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)

        gap_assay = []
        overlap_assay = []

        # sort table
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)
        
        # remove zero intervals 
        t=self.table[table_name][self.table[table_name]['TO']-self.table[table_name]['FROM']>=tol]


        if addgaps==True:
            # add gaps
            self.addtable(t,'__tmp',overwrite =True)
            gap_assay,overlap_assay = self.add_gaps(table_name='__tmp', 
                                                         new_table_name=new_table_name, 
                                                         overwrite = overwrite, 
                                                         tol=tol, 
                                                         endhole = endhole,
                                                         clean = clean)

            
            #remove table
            self.del_table('__tmp')
            
        else: 
            # just add table
            self.addtable(t,new_table_name,overwrite =True)
        
        # return outputs
        return gap_assay,overlap_assay
        

                

    # TODO: develop this: 
    # compositing
    
    cpdef collar2table(self, str table_name, 
                      str new_table_name,
                      collar_prop,
                      bint overwrite=False):
        """
                         
        mydhole.collar2table(table_name, new_table_name, 
                     collar_prop= ['TYPE', 'Comment'],
                     overwrite=False)                          
        
        Add collar properties to a table 
        
        
        Input:
        table_name: name of the table
        new_table_name: name of the output table
                
        collar_prop=: [list] with property names in collar 
            the property names may exist at drillhole.collar.columns
            if collar_prop= Null all the Collar properties will be copied
        
        overwrite: default False 
            If new_table_name exists and overwrite == True the existing 
            table will be overwrite.             
        
        """
        
        cdef int j, n
        
        assert table_name in self.table.keys(), 'The input table {} not in database'.format(table_name)
        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)
        
        if collar_prop == None: 
            
            out = pd.merge(self.table[table_name], self.collar,  on ='BHID')
        
        else: 
        
            assert isinstance (collar_prop, list), 'Error: collar_prop is not a list'
        
            n= len(collar_prop)
    
            # check that the properties exist in collar
            for p in collar_prop: 
                assert p in self.collar.columns, 'Error: The property {} not in collar'.format(p)
            
            p = ['BHID'] + collar_prop
            
            out = pd.merge(self.table[table_name], self.collar[p],  on ='BHID')
            
        
        self.addtable(table=out,table_name=new_table_name, overwrite=overwrite)
            
        
    
    cpdef downh_composite(self,  str table_name, 
                          str variable_name, 
                          str new_table_name, 
                          double cint = 1, 
                          double minlen=-1,                          
                          bint overwrite =False):

        """
        mydrillhole.downh_composite(
                          table_name, 
                          variable_name, 
                          new_table_name, 
                          cint = 1, 
                          minlen=-1,                          
                          overwrite =False)
        
        Down hole composite (one variable at the time)
        
               
        Input:
        table_name: name of the table with drillhole intervals
            it must be an existing table at drillhole.table.keys()
        variable_name: name of the variable in "table_name" table
            it must be an existing table at drillhole.table[table_name].columns
        new_table_name: name of the new table with composites
            it can be a new name, no in drillhole.table.keys()
            or a table name on drillhole.table.keys() if overwrite =True
        cint: default 1, compositing interval length
            it may be greater than zero 
        minlen: default -1 minimum composite length
            if < 0 then minlen will be cint/2.
        overwrite: default True. 
            If the table exist and overwrite = True the existing table 
            will be overwrite. 
        
        Notes
        ---
        Undefined intervals in 'variable_name' will be excluded. 
        
        To composite per lithology you can filter out each lithology
        composite and then combine the composite tables.
        
        This algorithm starts from distance zero and no residual will 
        be created at the top, for example, the if the first interval is 
        [FROM=1.5, TO=2] it will be composited as [FROM=0, TO=2,_len=0.5] 
        for 2 m composites. This produce consistent results if different
        tables are composited, making easier later table combination
        
        The accumulation output '_acum' can be use to produce accumulated
        variables, to identify number of intervals in each composite or 
        to track categoric variables (coded as number) compositing. 

        """
        
        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)
            
        #check that the variable is the table and the variable type
        assert variable_name in self.table[table_name].columns, 'Error: The variable {} do not exists in the input table'.format(variable_name)
        assert np.dtype(self.table[table_name][variable_name])==np.float64, 'Error: The variable {} is not type float64'.format(variable_name)
        
            
        
        # sort table
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)
                
        # create a group to easily iterate
        group=self.table[table_name].groupby('BHID')
        
        #add gaps
        BHID=group.groups.keys()
        nnf=[]
        nnt=[]
        nnBHID=[]
        nnvar= []
        nnacum = []
        nnlen = []
        for i in BHID:
                                                     
            nf, nt, nlen, nvar, nacum= composite1dh(ifrom= group.get_group(i)['FROM'].values,
                                                     ito= group.get_group(i)['TO'].values, 
                                                     ivar=group.get_group(i)[variable_name].values, 
                                                     cint=cint, minlen=minlen)
                                          
            
            nBHID = np.empty([len(nf)], dtype=object, order='C')
            nBHID[:]=i
            nnf+=nf.tolist()
            nnt+=nt.tolist()
            nnBHID+=nBHID.tolist()
            nnlen+=nlen.tolist()
            nnvar+=nvar.tolist()
            nnacum+=nacum.tolist()
    
        #create new table with gaps (only with fields )
        newtable=pd.DataFrame({'BHID':nnBHID, 'FROM':nnf,'TO':nnt,'_len':nnlen, variable_name:nnvar, '_acum':nnacum})
        
        #add table to the class
        self.addtable(newtable,new_table_name,overwrite)
                




    cpdef bench_composite(self, str table_name, 
                         double zmax,
                         double bench,
                         double tolerance=0.01):
        print 'we are working on that'        
    

    #-------------------------------------------------------------------
    #       VTK export 
    #-------------------------------------------------------------------
    cpdef export_core_vtk_line(self, str table_name, str filename,  double nanval=0, str title = '', bint binary=True):
        """
        export_core_vtk_line(table_name, filename, nanval=0, binary=True)
        
        Export desurveyed drillhole table to vtk lines. Endpoints 
        are required.
        
        
        Parameters
        ----------
        table_name : str
        filename : str 
            This is the absolute or relative file path and name. 
        title : str 
            The file header (or title, or comment)
        nanval: float
            Numeric fields with nan values will be replaced with nanval
        binary: boolean  (default True)
            If true export data in binary format, otherwise the output 
            will be in ascii format 
        
        
        See Also
        --------
        pygslib.vtktool 
        
        Notes
        -----
        Note: To export to VTK points see pygslib.vtktool 
        
        Examples
        --------
        
        >>> mydrillhole.export_core_vtk_line('assay', 'assay_line.vtk')
        
        """ 
        
        #check that table exists      
        assert table_name in self.table, '%s not exist in this drillhole database' % table_name
        
        #check columns 
        #check we have the rigth naming in collar columns 
        assert 'xe' in self.table[table_name].columns, "table without xe column"
        assert 'ye' in self.table[table_name].columns, "table without ye column"
        assert 'ze' in self.table[table_name].columns, "table without ze column"
        assert 'xb' in self.table[table_name].columns, "table without xb column"
        assert 'yb' in self.table[table_name].columns, "table without yb column"
        assert 'zb' in self.table[table_name].columns, "table without zb column"
        assert len(title)<=256, "The title exceeded the maximum 256 characters"

        cdef int l, n, dlen
        cdef np.ndarray[double, ndim=1] xb 
        cdef np.ndarray[double, ndim=1] yb 
        cdef np.ndarray[double, ndim=1] zb 
        cdef np.ndarray[double, ndim=1] xe 
        cdef np.ndarray[double, ndim=1] ye 
        cdef np.ndarray[double, ndim=1] ze 


        try: 
        
            # this will work if all coordinates defined
            
            assert  np.isfinite(self.table[table_name]['xb'].values).all(), "no finite coordinates at xb, please filter out non finite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['yb'].values).all(), "no finite coordinates at yb, please filter out non finite coordinates and try again" 
            assert  np.isfinite(self.table[table_name]['zb'].values).all(), "no finite coordinates at zb, please filter out non finite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['xe'].values).all(), "no finite coordinates at xe, please filter out non finite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['ye'].values).all(), "no finite coordinates at ye, please filter out non finite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['ze'].values).all(), "no finite coordinates at ze, please filter out non finite coordinates and try again"
            
            #create array views 
            xb = self.table[table_name]['xb'].values
            yb = self.table[table_name]['yb'].values
            zb = self.table[table_name]['zb'].values
            xe = self.table[table_name]['xe'].values
            ye = self.table[table_name]['ye'].values
            ze = self.table[table_name]['ze'].values        

        
        except:
            
            # this will work if there are coordinates non-defined (nan or non-finite)
            
            warnings.warn('There are non-finite or NAN values in coordinates, these intervals will be excluded')
            
            # using only xm as reference
            tmp = self.table[table_name][np.isfinite(self.table[table_name]['xm'])]
            tmp = tmp[np.isfinite(tmp['ym'])]
            tmp = tmp[np.isfinite(tmp['zm'])]
            #reset index (required to ensure right order???)
            tmp.reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')
            # a trick that makes this work, create a temporary drillhole
            self.addtable(tmp, '__tmp' ,overwrite = True)
            # update table name
            table_name = '__tmp'
            
            
            #create array views 
            xb = self.table[table_name]['xb'].values
            yb = self.table[table_name]['yb'].values
            zb = self.table[table_name]['zb'].values
            xe = self.table[table_name]['xe'].values
            ye = self.table[table_name]['ye'].values
            ze = self.table[table_name]['ze'].values   
 
            # if ther are undefined values at yb,zb, xe,ye,ze the output will be wrong, here whe check this
            assert  np.isfinite(self.table[table_name]['xb'].values).all(), "no finite coordinates at xb after nan at xm,ym,zm removed"
            assert  np.isfinite(self.table[table_name]['yb'].values).all(), "no finite coordinates at yb after nan at xm,ym,zm removed" 
            assert  np.isfinite(self.table[table_name]['zb'].values).all(), "no finite coordinates at zb after nan at xm,ym,zm removed"
            assert  np.isfinite(self.table[table_name]['xe'].values).all(), "no finite coordinates at xe after nan at xm,ym,zm removed"
            assert  np.isfinite(self.table[table_name]['ye'].values).all(), "no finite coordinates at ye after nan at xm,ym,zm removed"
            assert  np.isfinite(self.table[table_name]['ze'].values).all(), "no finite coordinates at ze after nan at xm,ym,zm removed"          
            
        
            
        
        #first we store the data in a set of vtk arrays (that will be cell data)
        dlen = xb.shape[0]
                
        vtkfields={}
        for i in self.table[table_name].columns:
            # assign the right vtk type 
            dtype = self.table[table_name][i].dtype
            if dtype==np.int: 
                vtkfields[i]= vtk.vtkIntArray()
            elif dtype==np.float or dtype==np.double: 
                vtkfields[i]= vtk.vtkDoubleArray()
            else:
                vtkfields[i]= vtk.vtkStringArray()
            
            # set some properties
            vtkfields[i].SetName(i)
            vtkfields[i].SetNumberOfComponents(1)
            vtkfields[i].SetNumberOfTuples(dlen)

            # deep copy data
            if dtype==np.int or dtype==np.float or dtype==np.double:
                for l in range(dlen): 
                    if np.isnan(self.table[table_name][i][l]):
                        vtkfields[i].SetValue(l,nanval)
                    else:
                        vtkfields[i].SetValue(l,self.table[table_name][i][l]) 
                        
            else:
                for l in range(dlen): 
                   vtkfields[i].SetValue(l,str(self.table[table_name][i][l]))
        
        # now we create a set of vtk points 
        points= vtk.vtkPoints()
        npoints = dlen*2

        # now we create a set of lines representing the cores and 
        # a line container (a cell array)
        line = vtk.vtkLine()
        lines = vtk.vtkCellArray()

        # populate this data
        n=-1
        for l in range(dlen):
            
            points.InsertNextPoint(xb[l], yb[l], zb[l])
            points.InsertNextPoint(xe[l], ye[l], ze[l])
            n=n+1
            line.GetPointIds().SetId(0,n)
            n=n+1
            line.GetPointIds().SetId(1,n)
            lines.InsertNextCell(line)
            
        
        # Create a polydata to store everything in
        linesPolyData = vtk.vtkPolyData()
         
        # Add the points to the dataset
        linesPolyData.SetPoints(points)
 
        # Add the lines to the dataset
        linesPolyData.SetLines(lines)
 
        #add properties 
        for i in vtkfields:
            linesPolyData.GetCellData().AddArray(vtkfields[i])

        # save data to VTK file
        writer = vtk.vtkPolyDataWriter()
        writer.SetInputData(linesPolyData)
        if binary == False: 
            writer.SetFileTypeToASCII()
        else: 
            writer.SetFileTypeToBinary()
        if title != '':
            writer.SetHeader (title)
        else: 
            writer.SetHeader ('Generated by PyGSLIB')
        writer.SetFileName(filename)
        writer.Write()
