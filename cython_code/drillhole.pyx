'''
PyGSLIB Drillhole, Module to handle drillhole data, desurvey
interval tables and other drillhole relate processes.

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
from libc.math cimport ceil
import pandas as pd
import warnings
import pygslib

#-------------------------------------------------------------------
#  General functions for desurvey
#-------------------------------------------------------------------
cpdef ang2cart( float azm,
                float dip):
    """ang2cart(float azm, float dip)

    Converts azimuth and dip to x, y, z.

    Returns the x,y,z coordinates of a 1 unit vector with origin of
    coordinates at 0,0,0 and direction defined by azm (azimuth) and
    dip (downward positive) angles.


    Parameters
    ----------
    azm,dip : float, in degrees
        The direction angles azimuth, with 0 or 360 pointing north and
        the dip angle measured from horizontal surface positive downward

    Returns
    -------
    out : tuple of floats, ``(x, y, z)``
        Cartesian coordinates.

    Examples
    --------
    >>> from pygslib.drillhole import ang2cart
    >>> ang2cart( azm = 45, dip = 75)
    (0.1830127239227295, 0.1830127090215683, -0.9659258127212524)
    >>>

    See Also
    --------
    cart2ang

    Note
    -----
    This function is to convert direction angle into a Cartesian
    representation, which are easy to interpolate. To convert
    back x,y,z values to direction angles use the function cart2ang.

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
    """cart2ang(float x, float y, float z)

    Converts x, y, z to azimuth and dip.

    Returns the azimuth and dip of a 1 unit vector with origin of
    coordinates at p1 [0,0,0] and p2 [x, y, z].


    Parameters
    ----------
    x,y,z : float
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

    Examples
    --------
    >>> from pygslib.drillhole import cart2ang
    >>> cart2ang(x = 0.18301, y = 0.18301, z = -0.96593)
    (45.0, 75.00092315673828)
    >>>

    See Also
    --------
    ang2cart


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

    azm= atan2(x,y)
    if azm<0.:
        azm= azm + pi*2
    azm = azm * RAD2DEG


    dip = -asin(z) * RAD2DEG

    return azm, dip


cpdef interp_ang1D( float azm1,
                    float dip1,
                    float azm2,
                    float dip2,
                    float len12,
                    float d1):
    """interp_ang1D( float azm1, float dip1, float azm2, float dip2, float len12, float d1)

    Interpolates the azimuth and dip angle over a line.

    Given a line with length ``len12`` and endpoints with direction
    angles ``azm1, dip1, azm2, dip2``, this function returns the
    direction angles of a point over this line, located at a distance
    d1 of the point with direction angles ``azm1, dip1``.


    Parameters
    ----------
    azm1,dip1,azm2,dip2,len12,d1 : float
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

    Example
    --------
    >>> from pygslib.drillhole import interp_ang1D
    >>> interp_ang1D(azm1=45, dip1=75, azm2=90, dip2=20, len12=10, d1=5)
    (80.74163055419922, 40.84182357788086)
    >>>

    See Also
    --------
    ang2cart, cart2ang

    Note
    -----
    The output direction angles are interpolated using average weighted
    by the distances ``d1`` and ``len12-d1``. To avoid issues with
    angles the azimuth and dip are converted to x,y,z, then interpolated
    and finally converted back to azimuth,dip angles.


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

cpdef __dsmincurb( float len12,
                 float azm1,
                 float dip1,
                 float azm2,
                 float dip2):

    """__dsmincurb( float len12, float azm1, float dip1, float azm2, float dip2)

    Desurveys one interval with minimum curvature

    Given a line with length ``len12`` and endpoints p1,p2 with
    direction angles ``azm1, dip1, azm2, dip2``, this function returns
    the differences in coordinate ``dz,dn,de`` of p2, assuming
    p1 with coordinates (0,0,0)

    Parameters
    ----------
    len12,azm1,dip1,azm2,dip2: float
        len12 is the length between a point 1 and a point 2.
        azm1, dip1, azm2, dip2 are direction angles azimuth, with 0 or
        360 pointing north and dip angles measured from horizontal
        surface positive downward. All these angles are in degrees.


    Returns
    -------
    out : tuple of floats, ``(dz,dn,de)``
        Differences in elevation, north coordinate (or y) and
        east coordinate (or x) in a Euclidean coordinate system.

    Example
    --------
    >>> from pygslib.drillhole import __dsmincurb
    >>> __dsmincurb(len12=10, azm1=45, dip1=75, azm2=90, dip2=20)
    (7.207193374633789, 1.0084573030471802, 6.186459064483643)
    >>>

    See Also
    --------
    ang2cart
    __dstangential

    Note
    -----
    The equations were derived from the paper

        `http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf`

    The minimum curvature is a weighted mean based on the
    dog-leg (dl) value and a Ratio Factor (rf = 2*tan(dl/2)/dl )
    if dl is zero we assign rf = 1, which is equivalent to  balanced
    tangential desurvey method. The dog-leg is zero if the direction
    angles at the endpoints of the desurvey intervals are equal.

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
        float rf
        float dl
        float DEG2RAD=3.141592654/180.0

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

    #print dz, dn, de, i1,a1,i2,a2,dl,rf

    return dz,dn,de


cpdef __dstangential( float len12,
                 float azm1,
                 float dip1):

    """__dstangential( float len12, float azm1, float dip1)

    Desurveys one interval with tangential method

    Given a line with length ``len12`` and endpoints p1,p2 with
    direction angle ``azm1, dip1``, this function returns
    the differences in coordinate ``dz,dn,de`` of p2, assuming
    p1 with coordinates (0,0,0)

    Parameters
    ----------
    len12,azm1,dip1: float
        len12 is the length between a point 1 and a point 2.
        azm1, dip1, are direction angle azimuth, with 0 or
        360 pointing north and dip angle measured from horizontal
        surface positive downward. All these angles are in degrees.


    Returns
    -------
    out : tuple of floats, ``(dz,dn,de)``
        Differences in elevation, north coordinate (or y) and
        east coordinate (or x) in a Euclidean coordinate system.

    Example
    -------
    >>> from pygslib.drillhole import __dstangential
    >>> __dstangential(len12=10, azm1=45, dip1=75)
    (9.659257888793945, 1.8301270008087158, 1.8301271200180054)
    >>>

    See Also
    --------
    ang2cart
    __dsmincurb


    Note
    -----
    The equations were derived from the paper

        `http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf`

    The tangential method has limitations but the error may be small
    for short intervals with direction angles interpolated from two
    desurvey points.


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
        float DEG2RAD = 3.141592654/180.0


    i1 = (90 - dip1) * DEG2RAD
    a1 = azm1 * DEG2RAD


    dz = len12*cos(i1)
    dn = len12*sin(i1)*cos(a1)
    de = len12*sin(i1)*sin(a1)

    return dz,dn,de

cpdef __angleson1dh(int indbs,
               int indes,
               np.ndarray[double, ndim=1] ats,
               np.ndarray[double, ndim=1] azs,
               np.ndarray[double, ndim=1] dips,
               float lpt,
               bint warns=True):
    """
    TODO: complete docstring
    """
    # output (angles at begin, mid and end interval)
    cdef:
        float azt
        float dipt

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

    assert ats[indbs]<EPSLON, 'first survey > 0 at %d' % indbs


    # for each point on desurvey
    for i in range (indbs,indes):



        # get the segment [a-b] to test interval
        a=ats[i]
        b=ats[i+1]
        azm1 = azs[i]
        dip1 = dips[i]
        azm2 = azs[i+1]
        dip2 = dips[i+1]
        len12 = ats[i+1]-ats[i]



        # test if we are in the interval, interpolate angles
        if lpt>=a and lpt<b:
            d1= lpt- a

            azt,dipt = interp_ang1D(azm1,dip1,azm2,dip2,len12,d1)

            return azt, dipt

    # no interva detected, we must be beyond the last survey!
    a=ats[indes]
    azt = azs[indes]
    dipt = dips[indes]
    # the point is beyond the last survey?
    if lpt>=a:

        return   azt, dipt

        if warns==True:
            warnings.warn('\n point beyond the last survey point at %s' % indes)
    else:
        if warns==True:
            warnings.warn('\n not interval found at survey, at %s' % indes)

        return   np.nan, np.nan


#-------------------------------------------------------------------
#  General functions for compositing
#-------------------------------------------------------------------
cpdef __composite1dh(double[:] ifrom,
                   double[:] ito,
                   double[:] ivar,
                   double cint = 1.,
                   double minlen=-1.):
    """__composite1dh(double[:] ifrom, double[:] ito, double[:] ivar, double cint = 1., double minlen=-1.)

    Composites intervals in a single drillhole. The From-To intervals may be sorted.

    Note
    ----
    See explanation of this implementation at
    http://opengeostat.com/downhole-compositing-algorithm/


    Parameters
    ----------
    ifrom,ito:     1D arrays of floats
        From - To  interval
    ivar :   1D array of floats
        variable to be composited
    cint: Optional, float, default 1.
        length of the compositing intervals
    minlen: Optional, float, default -1.
        minimum length of the composite, if <=0 then minlen = cint/2.

    Returns
    -------
    cfrom, cto:  1D arrays of floats
         From, To composited intervals
    clen, cvar, cacum:  1D arrays of floats
         total length of intervals composited
         variable composited
         variable accumulated

    Example
    -------
    >>>
    >>> cfrom, cto, clen, cvar, cacum= composite(ifrom,
                                                ito,
                                                ivar,
                                                cint = 1,
                                                minlen=-1)
    >>>

    See Also
    --------
    Drillhole.downh_composite

    TODO
    ----
     - [] Fix example section

    """


    assert len(ifrom) == len(ito)==len(ivar), "Error: ifrom, ito or ivar with different shapes"
    assert cint>0, "Error compositing length (cint) <= 0 "
    #assert all(ifrom < ito), "Error: ifrom >= ito, wrong or zero length intervals"
    #assert all(np.isfinite(ifrom)),"Error: ifrom with not finite elements"
    #assert all(np.isfinite(ito)),"Error: ito with not finite elements"

    cdef long i, l, k, ncomp, nintrb
    cdef double tmp


    #get some array parameters
    if minlen<0:
        minlen= cint/2.

    ncomp = int(ceil(ito[-1]/cint + 1))  # number of composites
    #ncomp = int(ito[-1]/cint + 1)  # number of composites
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
cdef __intersect(af,at,ai, bf,bt):
    ''' zf,zt,zi = __intersect(af,at,ai, bf,bt)

    Find intersect of intervals b on intervals a

    This is an auxiliar function to intersect drillhole tables
    
    Inputs
    ------
    af,at,ai: 1D, 1D numeric arrays, and 1D array 
        from, to, and id of intervals on table a
    bf,bt   : 1D, 1D numeric arrays 
        from, and to of intervals on table b
    Outputs
    -------
    zf,zt,zi arrays intersected

    Notes
    -----
    call this function with table a, b, and then b, a to alingn drillhole from to intervals

    Example
    -------
    >>>
    >>> ia, ib, l = __min_int(la, lb, ia, ib, tol=0.01)
    >>> af = [0,1,3,10,20,50] # from
    >>> at = [1,3,6,12,30,80] # to
    >>> ai = [0,1,2, 3, 4, 5] # id
    >>>
    >>> bf = [0,1,2,10]
    >>> bt = [1,2,3,70]
    >>> bi = [0,1,2, 3]
    >>> 
    >>> # intersect b on a
    >>> zf,zt,zi = intersect(af,at,ai, bf,bt)
    >>> b2a = pd.DataFrame({'FROM':zf, 'TO':zt, 'ID':zi})
    >>> 
    >>> # intersect a on b
    >>> zf,zt,zi = intersect(bf,bt,bi, af,at)
    >>> a2b = pd.DataFrame({'FROM':zf, 'TO':zt, 'ID':zi})
    >>>
    >>> # merge two tables
    >>> b2a.merge(a2b, on = ['FROM', 'TO'], how = 'outer').sort_values(by = ['FROM', 'TO'])

        FROM  TO  ID_x  ID_y
        0     0   1   0.0   0.0
        1     1   2   1.0   1.0
        2     2   3   1.0   2.0
        3     3   6   2.0   NaN
        4    10  12   3.0   3.0
        8    12  20   NaN   3.0
        5    20  30   4.0   3.0
        9    30  50   NaN   3.0
        6    50  70   5.0   3.0
        7    70  80   5.0   NaN

    >>>
    '''
    
    na = len(af)
    nb = len(bf)
    
    intersects = [set() for i in range(na)] # array of empty sets to contain intersects
    
    for i in range(na):
        for j in range(nb):
            # intersect FROM 
            if at[i]>bf[j] and af[i]<bf[j]:
                intersects[i].add(bf[j])
            # intersect TO 
            if at[i]>bt[j] and af[i]<bt[j]:
                intersects[i].add(bt[j])
    
    zf = []
    zt = []
    zi = []
    for i in range(na):
        zf.append(af[i])
        for j in sorted(intersects[i]):
            zt.append(j)
            zi.append(ai[i])
            zf.append(j)
        zt.append(at[i])
        zi.append(ai[i])
        
    return zf,zt,zi

cdef __min_int(double la,
             double lb,
             double ia,
             double ib,
             double tol=0.01):
    """__min_int(double la, double lb, double ia, double ib, double tol=0.01)

    This is an internal function used by __merge_one_dhole

    Given two complete drillholes A, B (no gaps and up to the end of
    the drillhole), this function returns the smaller of two
    intervals la = FromA[ia] lb = FromB[ib] and updates the
    indices ia and ib. There are three possible outcomes

    - FromA[ia] == FromB[ib]+/- tol. Returns mean of FromA[ia], FromB[ib] and ia+1, ib+1
    - FromA[ia] <  FromB[ib]. Returns FromA[ia] and ia+1, ib
    - FromA[ia] >  FromB[ib]. Returns FromB[ia] and ia, ib+1

    Example
    -------
    >>>
    >>> ia, ib, l = __min_int(la, lb, ia, ib, tol=0.01)
    >>>

    TODO
    ----
     - [] Fix example section

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

cdef __merge_one_dhole(double[:] la,
              double[:] lb,
              long[:] ida,
              long[:] idb,
              double tol=0.01):
    """ __merge_one_dhole(double[:] la, double[:] lb, long[:] ida, long[:] idb, double tol=0.01)

    Function to merge one drillhole.

    Note
    ----
    - Full drillhole is required (no gaps) and note that lengths are used
    instead From-To.
    - The length at the end of the drillhole must be the same

    Returns
    -------
    n: number of intervals
    lab: length from collar of each interval
    newida, newidb: Id linking data with the original tables a, b.

    Example
    -------

    >>>
    >>> n, lab, newida, newidb = __merge_one_dhole(la,lb, ida, idb, tol=0.01)
    >>>

    TODO
    ----
     - [] Fix example section
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
        ia, ib, l = __min_int(la[ia], lb[ib], ia, ib, tol=0.01)
        n+=1
        newida[n]=ida[ia-1]
        newidb[n]=idb[ib-1]
        lab[n]=l

        #this is the end of hole (this fails if maxdepth are not equal)
        if ia==maxia or ib==maxib:
            inhole=False

    return n, np_lab[:n+1], np_newida[:n+1], np_newidb[:n+1]

cdef __fillgap1Dhole(double[:] in_f,
            double[:] in_t,
            long[:] id,
            double tol=0.01,
            double endhole=-1):
    """__fillgap1Dhole(double[:] in_f, double[:] in_t, long[:] id, double tol=0.01, double endhole=-1)

    Function to fill gaps in one drillhole.

    Parameters
    ----------
    in_f,in_t,id: from, to intervals and interval ID.
        The interval ID is required to link back sample values after
        adding gaps
    tol: default 0.01. Tolerance
        gaps and overlaps <= tolerance will be ignored
    endhole: default -1. end of hold length
        if endhole>-1 a gap will be added if TO.last< endhole +/- tol
        if endhole>-1 and TO.last> endhole +/- tol a warning will be raised

    Returns
    -------
    np_nf,np_nt,np_nID: numpy 1D arrays with new from, to and interval ID
    np_gap,np_overlap: numpy 1D arrays with FROM position where gaps and
                       overlaps where detected.

    Example
    -------
    >>>
    >>> np_nf,np_nt,np_nID,np_gap,np_overlap = illgap(in_f, in_t, id, tol=0.01, endhole=-1)
    >>>

    TODO
    ----
     - [] Fix example section

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
                # this will happen only in unsorted array...
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



def groupcat(codes, table, ccol, cgroup, tablegrup= 'GROUP'):
    """

    Regroup categories using groping codes defined in a table.


    Parameters
    ----------
    codes: Pandas DataFrame
        DataFrame with categories and corresponding grouped categories
    table: Pandas DataFrame
        DataFrame with categories to be grouped
    ccol: str
        column with categories in DataFrame ``codes``
    cgroup: str
        column with grouped categories in DataFrame ``codes``
    tablegrup: str, defaul('GROUP')
        column with categories DataFrame ``table``


    Returns
    -------
    Dataframe with grouped categories
    Dict with grouping table

    TODO
    ----
     - [] Add example section

    """

    # make deep copy
    tmp = table.copy(deep=True)

    # get dictionary mapping types
    type_map = codes[[ccol, cgroup]].drop_duplicates().set_index(ccol)[cgroup].to_dict()

    # group categories
    tmp[tablegrup] = ""
    for i in type_map.keys():
        tmp.loc[tmp[ccol]== i, tablegrup] = type_map[i]

    return tmp, type_map



def txt2int(text):
    """txt2int(array text)

    Create an integer ID for cathegories in a text array

    Parameters
    ----------
    text : array of type str

    Returns
    ------
    array of integers

    TODO
    ----
     - [] Add example section

    """

    # create array of unique vales
    a= pd.unique(text) # this is a numpy array
    b= np.arange(len(a))
    c= pd.Series(b, index=a)

    return pd.Series(text).map(c)+1, c



#-------------------------------------------------------------------
#  Drillhole class
#-------------------------------------------------------------------
cdef class Drillhole:
    """Drillhole object with functions to desurvey and
    validate drillholes.


    Parameters
    ----------
    collar : Pandas DataFrame
    survey : Pandas DataFrame

    Example
    -------
    >>>
    >>> mydrillhole = pygslib.drillhole.Drillhole(collar, survey)
    >>> mydrillhole.addtable(assay, 'assay' ,overwrite = False)
    >>>


    The collar table may contain the fields:

     - BHID with any dtype
     - XCOLLAR, YCOLLAR, ZCOLLAR with dtypes `float64`.
     - (Optional) LENGTH with dtypes `float64`. This is the length of
       drillhole and can be used in some functions, for example,
       to fill gaps.

    The survey table may contain the fields:

      - BHID with any dtype
      - AT,AZ, DIP with dtypes `float64`.

    Attributes
    ----------
    collar : Pandas DataFrame
        The collar table
    survey : Pandas DataFrame
        The survey table
    tables : [Pandas DataFrame]
        list with Pandas DataFrames with interval tables (e.j. assay)
        containing compulsory fields: BHID with any dtype
        and FROM, TO with dtypes `float64`.
    table_mames : [str]
        list of table names

    Note
    ----
    A new copy of the input data will be created in memory. To work with
    shared memory in an external DataFrame you may copy back the table::

        >>> shared_collar = mydrillholeDB.collar

    - To add interval tables use ``addtable``
    - The existence of compulsory fields is only validated during
      object initialization and when adding interval tables
    - Survey tables may have at least 2 intervals and an interval AT=0
      is required for desurvey

    TODO
    ----
     - [] Fix example section.

    """
    cdef readonly object collar
    cdef readonly object survey
    cdef readonly object table

    property table_mames:
        def __get__(self):
            """
            Property getting
            """
            return self.table.keys()

    def __cinit__(self, collar, survey):
        """
        Class constructor

        Parameters
        ----------
        collar :  Pandas Dataframe

        survey :  Pandas Dataframe

        TODO
        ----
         - [] Add example section.

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

        # reset BHID to string
        self.collar['BHID'] = self.collar['BHID'].astype(str)
        self.survey['BHID'] = self.survey['BHID'].astype(str)

        # set BHID as uppercase to avoid issues
        self.collar['BHID']= self.collar['BHID'].str.upper()
        self.survey['BHID']= self.survey['BHID'].str.upper()

        # sort the data
        self.collar.sort_values(by=['BHID'], inplace=True)
        self.survey.sort_values(by=['BHID','AT'], inplace=True)


    cpdef addtable(self,object table,str table_name,bint overwrite =False):
        """addtable(object table,str table_name,bint overwrite =False)

        Adds a table and assigns a name.

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

        TODO
        ----
         - [] Fix example section.

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
        self.table[table_name]['BHID'] = self.table[table_name]['BHID'].astype(str)

        # set BHID as uppercase to avoid issues
        self.table[table_name]['BHID']= self.table[table_name]['BHID'].str.upper()

        # sort the data
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)

        # reset index, if necessary uncomment next line
        self.table[table_name].reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')

    cpdef del_table(self,str table_name):
        """del_table(str table_name)

        Deletes a table.

        Parameters
        ----------
        table_name : str
                the name of the table

        Examples
        --------

        >>> mydrillhole.addtable('assay')

        TODO
        ----
         - [] Fix example section.

        """
        #check the input is correct
        if table_name in self.table:
            del self.table[table_name]
        else:
            raise NameError('Table {} not found in the drillhole database object' % table_name)



    cpdef validate(self):
        """validate()

        Runs a set of basic validations on survey and collar
        consisting of:

        - check existence of Null values
        - check survey without values AT=0
        - check that coordinates and direction angles dtypes are float64
        - check survey without collar and collar without survey

        Returns
        -------
        dictionary with errors

        Examples
        --------
        >>>
        >>> mydrillhole.validate()
        >>>


        Note
        ----
        - You may run this validation before doing desurvey
        - Only few validation tests are implemented
        - Validation errors will raise a python error

        TODO
        ----
         - [] Implement check relation between table, large variations in
          survey angles, missing BHID.
         - [] Collect all errors and return a list of errors.
         - [] Collect all warnings and return a list of warnings.
         - [] Fix example section.
         - [] Check number range in survey (negative at and angles out of range)

        See Also
        --------
        validate_table

        """

        errors = {}

        #check collar

        #check repeated collar BHID
        if len(self.collar.loc[self.collar.duplicated(['BHID'])])>0:
            errors['collar duplicated'] = self.collar.loc[self.collar.duplicated(['BHID'])]['BHID'].values

        #check repeated survey BHID,AT
        if len(self.survey.loc[self.survey.duplicated(['BHID','AT'])])>0:
            errors['survey duplicated [BHID,AT]'] = self.survey.loc[self.survey.duplicated(['BHID','AT']), ['BHID','AT']].values

        # null values in collar
        if self.collar['BHID'].hasnans:
            errors['Non defined BHID in collar table'] = True
        if self.collar['XCOLLAR'].hasnans:
            errors['Non defined XCOLLAR in collar table'] = True
        if self.collar['YCOLLAR'].hasnans:
            errors['Non defined YCOLLAR in collar table'] = True
        if self.collar['ZCOLLAR'].hasnans:
            errors['Non defined ZCOLLAR in collar table'] = True
        if self.collar['XCOLLAR'].dtypes!='float64':
            errors['XCOLLAR in collar table != float64'] = True
        if self.collar['YCOLLAR'].dtypes!='float64':
            errors['YCOLLAR in collar table != float64'] = True
        if self.collar['ZCOLLAR'].dtypes!='float64':
            errors['ZCOLLAR in collar table != float64'] = True

        if 'LENGTH' in self.collar.columns:
            if self.collar['LENGTH'].dtypes!='float64':
                errors['LENGTH in collar table != float64'] = True
            if self.collar['LENGTH'].hasnans:
                errors['Non defined LENGTH values in collar table'] = True
        else:
            errors['No LENGTH column in collar'] = True


        #check SURVEY
        # null values in survey

        if self.survey['BHID'].hasnans:
            errors['Non defined BHID in survey table']= True
        if self.survey['AT'].hasnans:
            errors['Non defined AT in survey table']= True
        if self.survey['AZ'].hasnans:
            errors['Non defined AZ in survey table']= True
        if self.survey['DIP'].hasnans:
            errors['Non defined DIP in survey table']= True

        #check dtypes
        if self.survey['AT'].dtypes!='float64':
            errors['AT in survey table != float64']= True
        if self.survey['DIP'].dtypes!='float64':
            errors['DIP in survey table != float64']= True
        if self.survey['AZ'].dtypes!='float64':
            errors['AZ in survey table != float64']= True

        #Check using same data type for BHID
        if self.survey['BHID'].dtypes!=self.collar['BHID'].dtypes:
            errors["survey['BHID'].dtypes!=collar['BHID'].dtypes"]=True

        #check survey without values: at=0
        self.survey.sort_values(by=['BHID','AT'], inplace=True)
        error, bhid = self.__checkAt0(self.survey['BHID'].values, self.survey['AT'].values)
        if error>-1:
            errors['Firts interval AT!=0 at survey table, found at'] = bhid

        #check survey with only one survey occurrence per drillholes: this produces error in desurvey
        mask = self.survey.groupby(by='BHID').count()['AT'] == 1 # surveys with only one interval
        if mask.sum()>0:
          errors['Survey with one interval'] = self.survey.groupby(by='BHID').count()['AT'].loc[mask].index # drillholes with one interval

        #check survey without collar
        error = self.survey.loc[~self.survey['BHID'].isin(self.collar['BHID']), 'BHID'].values
        if error.shape[0]>0:
          errors['Survey without collar'] = error

        #check collar without survey
        error = self.collar.loc[~self.collar['BHID'].isin(self.survey['BHID']), 'BHID'].values
        if error.shape[0]>0:
          errors['Collar without survey'] = error

        # TODO: check survey.AT.last > endofhole

        return errors



    cdef __checkAt0(self,np.ndarray BHID, np.ndarray[double, ndim=1] AT):
        """
        This function checks if the first interval is approximately zero

        TODO
        ----
         - [] Fix example section.

        """
        # this is a hide function
        # the input data is assumed sorted
        cdef int n= AT.shape[0]
        cdef int i, start


        # first interval is zero
        if  AT[0]>0.00001:
            return 0, BHID[0]

        # the first dhole intervals (where BHID[i-1]!=BHID[i]) are zero?
        for i in range(1,n):
            if BHID[i-1]!=BHID[i]:
                if  AT[i]>0.00001:
                   return i, BHID[i]

        return -1, None


    cpdef fix_survey_one_interval_err(self, double dummy_at):
        """fix_survey_one_interval_err(self, double dummy_at)

        Fixes drillholes with a single survey record.

        Drillholes may have at least two survey records. This function
        adds a dummy survey record in drillholes with a single survey
        interval at a position ``dummy_at``.


        Parameters
        ----------
        dummy_at : float
            this value will be used as a dummy `AT` value at survey

        Example
        -------
        >>>
        >>> mydrillhole.fix_survey_one_interval_err(dummy_at = 900.0)
        >>>

        Note
        ----
        The dummy records added are a duplicated of the unique record
        available but with ``AT = dummy_at``. The result will be a
        straight drillhole.

        You will not be able to desurvey if there are drillholes with
        only one survey record.

        TODO
        ----
         - [] Fix example section.

        """

        mask = self.survey.groupby(by='BHID').count()['AT'] == 1 # surveys with only one interval
        surv_erro = self.survey.groupby(by='BHID').count()['AT'].loc[mask].index # drillholes with one interval
        survey_fix = self.survey.loc[self.survey['BHID'].isin(surv_erro)] # only those with one interval
        survey_fix['AT'] = dummy_at # dum at for desurvey
        self.survey = self.survey.append (survey_fix, ignore_index=True)

        #sort survey
        self.survey.sort_values(by=['BHID','AT'], inplace=True)


    cpdef validate_table(self,str table_name):
        """validate_table(str table_name)

        Runs a set of basic validations on interval tables
        consisting of:

        - checks null values in table BHID
        - checks null values in From/To
        - checks that FROM and TO dtypes are float64
        - checks relations between tables and collar

        Returns
        -------
        dictionary with errors

        See Also
        --------
        validate, add_gaps

        Examples
        --------

        >>> mydrillhole.validate_table('assay')


        """

        #check the input is correct
        assert table_name in self.table, '%s not exist in this drillhole database' % table_name

        #check table
        errors = {}

        #check repeated BHID,FROM
        if len(self.table[table_name][self.table[table_name].duplicated(['BHID','FROM'])])>0:
            tmp  =  self.table[table_name].loc[
                        self.table[table_name].duplicated(['BHID','FROM'], keep=False),
                        ['BHID','FROM']]
            errors['Duplicated BHID,FROM'] = tmp

        # null values in table bhid
        if self.table[table_name]['BHID'].hasnans:
            errors['Non defined BHID'] = True
        # null values in From/To
        if self.table[table_name]['FROM'].hasnans:
            errors['Non defined FROM'] = True
        if self.table[table_name]['TO'].hasnans:
            errors['Non defined TO'] = True
        if self.table[table_name]['FROM'].dtypes!='float64':
            errors['FROM in table != float64' ] = True
        if self.table[table_name]['TO'].dtypes!='float64':
            errors['TO in table != float64' ] = True

        #Check using same data type for BHID
        if self.table[table_name]['BHID'].dtypes!=self.collar['BHID'].dtypes:
            errors["Table['BHID'].dtypes!=Collar['BHID'].dtypes"] = True

        #check table without collar
        tmp = self.table[table_name].loc[~self.table[table_name]['BHID'].isin(self.collar['BHID'].unique()), 'BHID'].values
        if tmp.shape[0]>0:
          errors['Table BHID not at collar BHID list'] = tmp

        #check collar without table
        tmp = self.collar.loc[~self.collar['BHID'].isin(self.table[table_name]['BHID'].unique()), 'BHID'].values
        if tmp.shape[0]>0:
          errors['Collar BHID not at table BHID list'] = tmp


        # TODO: check overlaps. Done see add gaps
        # TODO: check TO.last > endofhole

        return errors

    cpdef txt2intID(self, str table_name):
        """txt2intID(str table_name)

        Creates an alternative BHID of type integer

        A new ``BHIDint`` will be created on the ``table[table_name]``
        and in collar. ``BHIDint`` is just and ordered list of integers.

        BHIDint may be required in some functions compiles in Fortran,
        for example ``pyslib.gslib.gamv`` and ``pyslib.gslib.kt3d``.

        Parameters
        ----------
        table_name : str

        Examples
        --------
        >>>
        >>> mydrillhole.collar.sort(['BHID'], inplace=True)
        >>> mydrillhole.table['assay'].sort(['BHID', 'FROM'], inplace=True)
        >>> mydrillhole.txt2intID('assay')
        >>>

        Note
        ----
        The Collar and the table may be sorted.

        TODO
        ----
         - [] Fix example section.

        """

        # the data is assumed sorted

        assert table_name in self.table, 'The table %s do not exist in the database' % table_name

        cdef int sloc
        cdef int i
        cdef int j
        cdef int nc = self.collar.shape[0]
        cdef int nt = self.table[table_name].shape[0]
        cdef np.ndarray[long, ndim=1] cBHID = np.zeros([nc], dtype=int)
        cdef np.ndarray[long, ndim=1] tBHID = np.zeros([nt], dtype=int)
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


    cpdef desurvey_survey(self, method = 1):
        """desurvey_survey()
        Add coordinates to the survey table. 
        
        Parameters
        ----------
        method : int, optional, default 1 (minimum curvature)
            the desurvey method: 1 is minimum curvature any other value is tangential
        Note
        -----
        Coordinates in survey table are used to desurvey interval tables. 

        This function will be call from desurvey2 if the fields x,y, 
        and z ar not in the survey table, using the same `method`. 

        Make sure to use the same desurvey technique to 
        desurvey interval and survey tables.

        """
        self.survey['x'] = self.survey['y'] = self.survey['z'] = np.nan

        self.survey.sort_values(by = ['BHID', 'AT'])

        for c in self.collar['BHID']:
            # find collar data for this drillhole
            XC,YC,ZC = self.collar.loc[self.collar['BHID']==c, ['XCOLLAR','YCOLLAR','ZCOLLAR']].values[0]

            # desurvey survey
            AT = self.survey.loc[self.survey['BHID']==c, 'AT'].values
            DIP = self.survey.loc[self.survey['BHID']==c, 'DIP'].values
            AZ = self.survey.loc[self.survey['BHID']==c, 'AZ'].values

            dz = np.empty(AT.shape)
            dn = np.empty(AT.shape)
            de = np.empty(AT.shape)

            x = np.empty(AT.shape)
            y = np.empty(AT.shape)
            z = np.empty(AT.shape)

            dz[0] = dn[0] = de[0] = 0
            x[0] = XC
            y[0] = YC
            z[0] = ZC
            
            
            for i in range(1, AT.shape[0]):
                if method == 1:
                    dz[i],dn[i],de[i] = __dsmincurb(len12 = AT[i] - AT[i-1], azm1 = AZ[i-1],  dip1 = DIP[i-1], azm2=AZ[i], dip2 = DIP[i])
                else:
                    dz[i],dn[i],de[i] = __dstangential(len12 = AT[i] - AT[i-1], azm1 = AZ[i-1],  dip1 = DIP[i-1])

                x[i] = x[i-1] + de[i] 
                y[i] = y[i-1] + dn[i]
                z[i] = z[i-1] - dz[i]
                
            self.survey.loc[self.survey['BHID']==c,'x'] = x
            self.survey.loc[self.survey['BHID']==c,'y'] = y
            self.survey.loc[self.survey['BHID']==c,'z'] = z     

    cpdef desurvey_table(self, str table_name, int method=1):
        """desurvey_table(str table_name, int method=1)

        Desurvey a drillhole table.

        Create coordinates and direction angles at table midpoint
        intervals. If ``endpoints=True`` it also creates coordinate fields
        at end point intervals. The existing coordinate fields will be
        overwritten.

        Parameters
        ----------
        table_name : str
            a table name existing in drillhole object
        method : int, optional, default 1 (minimum curvature)
            the desurvey method: 1 is minimum curvature and 2 is tangential


        Examples
        --------
        >>> mydrillhole.desurvey('assay', method=1)
        >>>

        Note
        ----
        This function calls __dsmincurb() or
        __dstangential() functions to calculate the desurvey value.
        
        If the last survey interval is shallower than the deepest point in 
        the table, then the las survey interval is reapeated at the depth of
        the deepest point in the table + 0.001. 

        Both desurvey methods (tangential and minimum curvature) use
        angles interpolated from two desurvey points at Survey table.

        TODO
        ----
         - [] Fix example section.
         - [] Optimize

        """

        # first desurvey survey
        if 'x' not in self.survey.columns or \
           'y' not in self.survey.columns or \
           'z' not in self.survey.columns:

           self.desurvey_survey(method = method)          

        # prepare output (this is a deep copy of the table (?))
        table = self.table[table_name].set_index('BHID')
        table ['xb'] = np.nan
        table ['yb'] = np.nan
        table ['zb'] = np.nan
        table ['xe'] = np.nan
        table ['ye'] = np.nan
        table ['ze'] = np.nan
        table ['xm'] = np.nan
        table ['ym'] = np.nan
        table ['zm'] = np.nan
        table ['azmb'] = np.nan
        table ['dipb'] = np.nan
        table ['azme'] = np.nan
        table ['dipe'] = np.nan
        table ['azmm'] = np.nan
        table ['dipm'] = np.nan

        survey = self.survey.set_index('BHID')

        for c in self.table[table_name]['BHID'].unique(): 
            
            # get survey
            AT =  survey.loc[c, ['AT']].values.ravel()  
            DIP = survey.loc[c, ['DIP']].values.ravel()
            AZ =  survey.loc[c, ['AZ']].values.ravel()
            xs =  survey.loc[c, ['x']].values.ravel()
            ys =  survey.loc[c, ['y']].values.ravel()
            zs =  survey.loc[c, ['z']].values.ravel()

            # get from, to, y mid interval
            db = table.loc[c, ['FROM']].values.ravel() # [[]].values.ravel() is to prevent getting a scalar if shape is 1
            de = table.loc[c, ['TO']].values.ravel()
            dm = db + (de-db)/2

            # add at the end of the survey if de< AT
            if de[-1]>= AT[-1]:
                AZ = np.append(AZ, AZ[-1])
                DIP = np.append(DIP, DIP[-1])
                AT = np.append(AT, de[-1] + 0.01)

            #get the index where each interval is located
            jb = np.searchsorted(AT, db, side='right')
            je = np.searchsorted(AT, de, side='right')
            jm = np.searchsorted(AT, dm, side='right')

            # outputs
            azmt = np.empty(jb.shape)
            dipt = np.empty(jb.shape)           
            x = np.empty(jb.shape)
            y = np.empty(jb.shape)
            z = np.empty(jb.shape)

            # the bigining
            for i in range(jb.shape[0]):
                d1 = db[i] -AT[jb[i]-1]
                lll1 = AT[jb[i]]
                lll2 = AT[jb[i]-1]
                len12 = lll1-lll2
                azm1 = AZ[jb[i]-1]
                dip1 = DIP[jb[i]-1]
                azm2 = AZ[jb[i]]
                dip2 = DIP[jb[i]]
                azmt[i],dipt[i] = interp_ang1D(azm1, dip1, azm2, dip2, len12, d1)
                if method==1:
                    dz,dy,dx = __dsmincurb(d1, azm1,  dip1, azmt[i], dipt[i])
                else:
                    dz,dy,dx = __dstangential(d1, azm1,  dip1)
                
                x[i] = dx + xs[jb[i]-1]
                y[i] = dy + ys[jb[i]-1]
                z[i] = zs[jb[i]-1] - dz
    
            table.loc[c,'azmb']  = azmt
            table.loc[c,'dipb']  = dipt
            table.loc[c,'xb']  = x
            table.loc[c,'yb']  = y
            table.loc[c,'zb']  = z

            # the end
            for i in range(je.shape[0]):
                d1 = de[i] -AT[je[i]-1]
                len12 = AT[je[i]]-AT[je[i]-1]
                azm1 = AZ[je[i]-1]
                dip1 = DIP[je[i]-1]
                azm2 = AZ[je[i]]
                dip2 = DIP[je[i]]
                azmt[i],dipt[i] = interp_ang1D(azm1, dip1, azm2, dip2, len12, d1)
                if method==1:
                    dz,dy,dx = __dsmincurb(d1, azm1,  dip1, azmt[i], dipt[i])
                else:
                    dz,dy,dx = __dstangential(d1, azm1,  dip1)
                x[i] = dx + xs[je[i]-1]
                y[i] = dy + ys[je[i]-1]
                z[i] = zs[je[i]-1] - dz
                 
            table.loc[c,'azme']  = azmt
            table.loc[c,'dipe']  = dipt
            table.loc[c,'xe']  = x
            table.loc[c,'ye']  = y
            table.loc[c,'ze']  = z

            # the mean 
            for i in range(jm.shape[0]):
                d1 = dm[i] -AT[jm[i]-1]
                len12 = AT[jm[i]]-AT[jm[i]-1]
                azm1 = AZ[jm[i]-1]
                dip1 = DIP[jm[i]-1]
                azm2 = AZ[jm[i]]
                dip2 = DIP[jm[i]]
                azmt[i],dipt[i] = interp_ang1D(azm1, dip1, azm2, dip2, len12, d1)
                if method==1:
                    dz,dy,dx = __dsmincurb(d1, azm1,  dip1, azmt[i], dipt[i])
                else:
                    dz,dy,dx = __dstangential(d1, azm1,  dip1)

                x[i] = dx + xs[jm[i]-1]
                y[i] = dy + ys[jm[i]-1]
                z[i] = zs[jm[i]-1] - dz
                
               
            table.loc[c,'azme']  = azmt
            table.loc[c,'dipm']  = dipt
            table.loc[c,'xm']  = x
            table.loc[c,'ym']  = y
            table.loc[c,'zm']  = z
          
        # update
        self.table[table_name] = table.reset_index()

        # this is to produce a warning if some intervals where not desurvey
        try:
            assert  np.isfinite(self.table[table_name]['xm'].values).all()
            assert  np.isfinite(self.table[table_name]['xb'].values).all()
            assert  np.isfinite(self.table[table_name]['xe'].values).all()
        except:
            warnings.warn('Some intervals where non-desurveyed and NAN coordinates where created, check errors in collar or survey')

    cpdef desurvey(self, str table_name, bint endpoints=True,
                   bint warns=True, int method=1):
        """desurvey(str table_name, bint endpoints=True, bint warns=True, int method=1)

        *DEPRECATED* 

        Desurvey a drillhole table.

        Create coordinates and direction angles at table midpoint
        intervals. If ``endpoints=True`` it also creates coordinate fields
        at end point intervals. The existing coordinate fields will be
        overwritten.

        .. deprecated:: 
            Use desrvey table instead, this fuction has bugs. It produce a "small" deviatio that can be 
            important in some special cases. 

        Parameters
        ----------
        table_name : str
            a table name existing in drillhole object
        endpoints : boolean, optional, default True
            if True coordinates at beginning and end intervals are also created
        warns : boolean, optional, default True
            if True warns will be raised if inconsistencies are detected
        method : int, optional, default 1 (minimum curvature)
            the desurvey method: 1 is minimum curvature and 2 is tangential


        Examples
        --------
        >>> # sort with pandas function DataFrame.sort
        >>> mydrillhole.collar.sort(['BHID'], inplace=True)
        >>> mydrillhole.survey.sort(['BHID', 'AT'], inplace=True)
        >>> mydrillhole.table['assay'].sort(['BHID', 'FROM'], inplace=True)
        >>> # desurvey
        >>> mydrillhole.desurvey('assay', endpoints=True, method=1)
        >>>


        Note
        ----
        `collar`, `survey` and the input table may be sorted.
        If you call the function with `endpoints=False` end points already
        desurveyed may not be overwritten.

        This function calls __dsmincurb() and __desurv1dh() or
        __dstangential() functions to calculate the desurvey value.

        Both desurvey methods (tangential and minimum curvature) use
        angles interpolated from two desurvey points at Survey table.

        TODO
        ----
         - [] Fix example section.

        """

        warnings.warn('This function has bugs and is deprecated, use desurvey_table instead')


        # check the input is correct
        assert table_name in self.table, "table %s not exist" % table_name
        assert method==1 or method==2, " invalid parameter method=%, valid inputs are 1 or 2 " % method


        # sort tables
        self.collar.sort_values(by=['BHID'], inplace=True)
        self.survey.sort_values(by=['BHID', 'AT'], inplace=True)
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)


        # define arrays
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
        cdef double mid
        cdef double de
        cdef double dn
        cdef double dz
        cdef double x
        cdef double y
        cdef double z
        cdef double azm1
        cdef double azm2
        cdef double dip1
        cdef double dip2
        cdef double at

        # otput
        cdef np.ndarray[double, ndim=1] azmt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] dipmt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] xmt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] ymt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] zmt = np.empty([nt], dtype=float)


        #if endpoints==true:
        cdef np.ndarray[double, ndim=1] azbt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] dipbt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] xbt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] ybt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] zbt = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] azet = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] dipet = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] xet = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] yet = np.empty([nt], dtype=float)
        cdef np.ndarray[double, ndim=1] zet = np.empty([nt], dtype=float)

        # initialize output to nans
        azmt[:] = np.nan
        dipmt[:] = np.nan
        azbt[:] = np.nan
        dipbt[:] = np.nan
        azet[:] = np.nan
        dipet[:] = np.nan
        xmt[:] = np.nan
        ymt[:] = np.nan
        zmt[:] = np.nan
        xbt[:] = np.nan
        ybt[:] = np.nan
        zbt[:] = np.nan
        xet[:] = np.nan
        yet[:] = np.nan
        zet[:] = np.nan


        # first pass calculate angles and de,dn,dz

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

            # if survey not found then skip this drillhole
            if indbs==-1 or indes==-1:
                # do not desurvey this drillhole
                warnings.warn('! collar {} without survey, table not desurveyed'.format(idc[jc]))
                continue


            # angles over 1st desurvey point == collar
            azm1  = azs[indbs]
            dip1 = dips[indbs]
            at = 0.


            # TODO: remove this and allow desurvey with only one record

            # if only one survey points then skip this dillhole
            if indbs==indes:
                # do not desurvey this drillhole
                warnings.warn('! collar {} without survey at end collar, table not desurveyed'.format(idc[jc]))
                continue

            # initialize coordinates
            x =  xc[jc]
            y =  yc[jc]
            z =  zc[jc]

            # for each point in the table
            for jt in range(indt, nt):
                # the table id is similar to collar? Then desurvey


                if idc[jc]==idt[jt]:

                    indt = jt # do not loop again before this index

                    # a) beginning interval
                    azm2,dip2 = __angleson1dh(indbs,indes,ats,azs,dips,fromt[jt],warns)
                    azbt[jt] = azm2
                    dipbt[jt] = dip2
                    len12 = float(fromt[jt]) - at

                    if method == 1 :
                        dz,dn,de = __dsmincurb(len12,azm1,dip1,azm2,dip2)
                    elif method == 2 :
                        dz,dn,de = __dstangential(len12,azm1,dip1)

                    xbt[jt] = de
                    ybt[jt] = dn
                    zbt[jt] = dz


                    # update for next interval
                    azm1 = azm2
                    dip1 = dip2
                    at   = float(fromt[jt])


                    # b) mid interv

                    mid = float(fromt[jt]) + float(tot[jt]-fromt[jt])/2.

                    azm2, dip2 = __angleson1dh(indbs,indes,ats,azs,dips,mid,warns)
                    azmt[jt] = azm2
                    dipmt[jt]= dip2
                    len12 = mid - at

                    if method == 1 :
                        dz,dn,de = __dsmincurb(len12,azm1,dip1,azm2,dip2)
                    elif method == 2 :
                        dz,dn,de = __dstangential(len12,azm1,dip1)

                    xmt[jt] = de + xbt[jt]
                    ymt[jt] = dn + ybt[jt]
                    zmt[jt] = dz + zbt[jt]

                    # update for next interval
                    azm1 = azm2
                    dip1 = dip2
                    at   = mid

                    # c) end interval

                    azm2, dip2 = __angleson1dh(indbs,indes,ats,azs,dips,float(tot[jt]),warns)
                    azet[jt] = azm2
                    dipet[jt] = dip2
                    len12 = float(tot[jt]) - at

                    if method == 1 :
                        dz,dn,de = __dsmincurb(len12,azm1,dip1,azm2,dip2)
                    elif method == 2 :
                        dz,dn,de = __dstangential(len12,azm1,dip1)

                    xet[jt] = de + xmt[jt]
                    yet[jt] = dn + ymt[jt]
                    zet[jt] = dz + zmt[jt]

                    # update for next interval
                    azm1 = azm2
                    dip1 = dip2
                    at   = float(tot[jt])


                    # now we calculate coordinates
                    xbt[jt] = x+float(xbt[jt])
                    ybt[jt] = y+float(ybt[jt])
                    zbt[jt] = z-float(zbt[jt])
                    xmt[jt] = x+float(xmt[jt])
                    ymt[jt] = y+float(ymt[jt])
                    zmt[jt] = z-float(zmt[jt])
                    xet[jt] = x+float(xet[jt])
                    yet[jt] = y+float(yet[jt])
                    zet[jt] = z-float(zet[jt])

                    #print jt, xbt[jt],ybt[jt],zbt[jt]
                    #print jt, xmt[jt],ymt[jt],zmt[jt]
                    #print jt, xet[jt],yet[jt],zet[jt]

                    # update for next interval
                    x = xet[jt]
                    y = yet[jt]
                    z = zet[jt]


        self.table[table_name]['azmm'] = azmt
        self.table[table_name]['dipm']= dipmt
        self.table[table_name]['xm']= xmt
        self.table[table_name]['ym']= ymt
        self.table[table_name]['zm']= zmt
        if endpoints==True:
            self.table[table_name]['azmb'] = azbt
            self.table[table_name]['dipb']= dipbt
            self.table[table_name]['xb']= xbt
            self.table[table_name]['yb']= ybt
            self.table[table_name]['zb']= zbt
            self.table[table_name]['azme'] = azet
            self.table[table_name]['dipe']= dipet
            self.table[table_name]['xe']= xet
            self.table[table_name]['ye']= yet
            self.table[table_name]['ze']= zet



        # this is to produce a warning if some intervals where not desurvey
        try:
            assert  np.isfinite(self.table[table_name]['xm'].values).all(), "nofinite coordinates at xb, please filter out nonfinite coordinates and try again"
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

        """add_gaps(str table_name, str new_table_name, bint overwrite =False, double tol=0.01, bint endhole=False, bint clean=True)

        Fills gaps with new FROM-TO intervals.

        All the gaps, including the gaps at collar and at the end of
        drillholes will be filled with ``NaN`` intervals.

        A code field ``_id0`` will be added to new and existing tables.
        ``_id0`` can be used to link input table rows with the output
        table rows. Gaps will have ``_id0= -999``


        Parameters
        ----------
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

        Returns
        -------
        gap,overlap: list


        Example
        -------
        >>>
        >>> gap,overlap= mydrillhole.add_gaps(table_name = 'assay',
                                             new_table_name = 'tmp',
                                             overwrite=False,
                                             tol=0.01,
                                             endhole=False,
                                             clean=True)
        >>>


        Note
        ----
        The output ``gap`` and ``overlap`` contain ``_id0`` of the raws
        before the gap or the overlap was detected.
        Gaps and overlaps at the end of the drillhole
        (if endhole==True) _id0 will have value -888.

        This function only works if there is more than one drillhole in collar

        TODO
        ----
         - [] Fix example section.


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

            nf,nt,nID,gap,overlap=__fillgap1Dhole(in_f = group.get_group(i)['FROM'].values,
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


    cpdef split_long_intervals(self,
                     str table,
                     str new_table,
                     double maxlength,
                     double splitlength=0,
                     double minlength=0.,
                     bint overwrite =False,
                     bint clean=True,
                     long buff=10):
        """split_long_intervals(str table, str new_table, double maxlength, double splitlength, double minlength=0, bint overwrite =False, bint clean=True, long buff=10)

        Split long sample intervals in a table.

        This function split sample intervals longer than ``maxlength`` into
        intervals approximatelly ``splitlength``, to be precise, equal to

          ``newlength = interval_lenght / splitlength``

        The interval will be not splitted if ``newlength < minlength``



        Parameters
        ----------
        table: name of the first table
            it must be an existing table at drillhole.table.keys()
        new_table: name of the new table
            it may not exists at drillhole.table.keys() if overwrite == False
        maxlength:
            all intervals longer than this length will be splitted if
            the resulting splits are longuer than minlength
        splitlength: default 0.
            reference split length. If <= 0 or not specified it will be
            set equal to maxlength
        minlength: default 0.
            intervals will not be splitted if the resulting split is
            less than this value. If <= 0 or not specified it will be
            set equal to splitlength/2
        overwrite: default True.
            If new_table_name exists and overwrite == True the existing
            table will be overwritten.
        clean: default True.
            Delete temporary columns created with suffix __tmp__.
        buff: defaul 10
            an internal and temporary table with 10*total_length/splitlength
            interval will be created in memory. If you have an error
            in memory use a larger number.


        Example
        -------

        >>> mydrillhole.split_long_intervals(
                     table= "litho",
                     new_table = "litho_split",
                     maxlength = 5,
                     splitlength=4.5,
                     minlength=2.,
                     overwrite =False,
                     clean=True)
        >>>
        >>>


        TODO
        ----
         - [] Fix example section.

        """



        # check that the table is not in the database
        if overwrite==False:
            assert new_table not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table)

        # check splitlength
        if splitlength < 0.000001:
            splitlength = maxlength

        # check minlength
        if minlength < 0.000001:
            minlength = splitlength / 2.

        # check
        assert splitlength <= maxlength, 'Error in input parameters splitlength > maxlength'
        assert minlength <= splitlength, 'Error in input parameters minlength > splitlength'



        # add ID to table and sort
        self.table[table].loc[:,'_id0']= np.arange(self.table[table].shape[0], dtype=long)[:]
        self.table[table].loc[:,'_LENGTH']= self.table[table]['TO']-self.table[table]['FROM']
        self.table[table].sort_values(by=['BHID', 'FROM'], inplace=True)


        # define arrays
        cdef np.ndarray[long, ndim=1] idt =self.table[table]['_id0'].values
        cdef np.ndarray[double, ndim=1] fromt = self.table[table]['FROM'].values
        cdef np.ndarray[double, ndim=1] tot = self.table[table]['TO'].values
        cdef np.ndarray[double, ndim=1] length = self.table[table]['_LENGTH'].values


        cdef:
             double splitl
             long ll
             long n
             int nsplit

        ll = length.sum()/ splitlength * buff


        new_id0 = np.empty(ll, dtype=np.int64, order='C')
        new_from = np.empty(ll, dtype=np.float64, order='C')
        new_to = np.empty(ll, dtype=np.float64, order='C')


        n = -1
        for i in range(self.table[table].shape[0]):
            # split if >= maxleng
            if length[i] >= maxlength:
                # the new split is ok
                nsplit = <int> (length[i] / splitlength + 0.5)    # this may be round properly
                splitl = length[i]  / nsplit


                if splitl >= minlength:
                    for j in range(nsplit):
                        # split here
                        n+=1
                        new_id0[n]= idt[i]
                        new_from[n]= fromt[i] + j*splitl
                        new_to[n]= new_from[n] + splitl

                else:
                     n+=1
                     if n<ll:
                         new_id0[n]= idt[i]
                         new_from[n]= fromt[i]
                         new_to[n]= tot[i]

            else:
                n+=1
                if n<ll:
                    new_id0[n]= idt[i]
                    new_from[n]= fromt[i]
                    new_to[n]= tot[i]
                else:
                    # classsic python error
                    assert n<ll, 'Error allocating internal table in memory. Use a larger buff value'


        n+=1


        #create new table with intervals and ID
        newtable=pd.DataFrame({'FROM':new_from[:n],'TO':new_to[:n],'_id0':new_id0[:n]})


        # merge with existing data and sort
        newtable=newtable.join(self.table[table], on='_id0', rsuffix='__tmp__')
        newtable.sort_values(by=['BHID', 'FROM'], inplace=True)

        #clean if necessary
        if clean:
            newtable.drop(
               ['FROM__tmp__','TO__tmp__','_id0__tmp__'],
               axis=1,inplace=True, errors='ignore')

        #add table to the class
        self.addtable(newtable,new_table,overwrite)




    cpdef merge(self,str table_A,
                     str table_B,
                     str new_table_name,
                     bint overwrite =False,
                     double tol=0.01,
                     bint clean=True):
        """Deprecated, this function has bugs use merge_table() instead
        merge(str table_A, str table_B, str new_table_name, bint overwrite =False, double tol=0.01, bint clean=True)

        Combines two tables by intersecting intervals.

        This function requires drillholes without gaps and overlaps.
        You may un add_gaps in table_A and table_B before using
        this function.


        Parameters
        ----------
        table_A: name of the first table
            it must be an existing table at drillhole.table.keys()
        table_B: name of the second table
            it must be an existing table at drillhole.table.keys()
        new_table_name: name of the new table
            it may not exists at drillhole.table.keys()
        overwrite: default True.
            If new_table_name exists and overwrite == True the existing
            table will be overwritten.
        tol: default 0.01.
            segments, gaps and overlays within tol will be ignore but
            adjusted
        clean: default True.
            Delete temporary columns created with suffix __tmp__.


        Example
        -------

        >>> mydrillhole.merge(table_A = 'assay',
                            table_B = 'litho',
                            new_table_name = 'tmp',
                            overwrite =False,
                            tol=0.01,
                            clean=True)
        >>>
        >>>

        TODO
        ----
         - [] Fix example section.

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
        np_la=np.empty(self.table[table_A].shape[0]+1, dtype = float)
        np_lb=np.empty(self.table[table_B].shape[0]+1, dtype = float)
        np_ida=np.empty(self.table[table_A].shape[0]+1, dtype = int)
        np_idb=np.empty(self.table[table_B].shape[0]+1, dtype = int)

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
                n, np_lab, np_newida, np_newidb = __merge_one_dhole(la[:nk+1],lb[:nj+1], ida[:nk+1], idb[:nj+1], tol=0.01)

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


    cpdef merge_table(self,str table_A,
                     str table_B,
                     str new_table_name,
                     bint overwrite =False,
                     bint clean=True):
        """merge_table(str table_A, str table_B, str new_table_name, bint overwrite =False, bint clean=True)

        Combines two tables by intersecting intervals.


        Parameters
        ----------
        table_A: name of the first table
            it must be an existing table at drillhole.table.keys()
        table_B: name of the second table
            it must be an existing table at drillhole.table.keys()
        new_table_name: name of the new table
            it may not exists at drillhole.table.keys()
        overwrite: default True.
            If new_table_name exists and overwrite == True the existing
            table will be overwritten.
        clean: default True.
            Delete temporary columns created with suffix __tmp__.


        Example
        -------

        >>> mydrillhole.merge(table_A = 'assay',
                            table_B = 'litho',
                            new_table_name = 'tmp',
                            overwrite =False,
                            clean=True)
        >>>
        """


        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)


        # sort tables
        self.table[table_A].sort_values(by=['BHID', 'FROM'], inplace=True)
        self.table[table_B].sort_values(by=['BHID', 'FROM'], inplace=True)

        # add ID to tables
        self.table[table_A].loc[:,'_ida']= np.arange(self.table[table_A].shape[0])[:]
        self.table[table_B].loc[:,'_idb']= np.arange(self.table[table_B].shape[0])[:]

        # create a groups to easily iterate
        groupA=self.table[table_A].groupby('BHID')
        groupB=self.table[table_B].groupby('BHID')

        keysA= groupA.groups.keys()
        keysB= groupB.groups.keys()


        #merge
        BHID=self.collar.BHID.values
        taf=[]
        tat=[]
        tai=[]
        ta_bhid=[]
        tbf=[]
        tbt=[]
        tbi=[]
        tb_bhid=[]

        for i in BHID:
            if (i not in keysA) and (i not in keysB):
                continue
            #get arrays
            if i in self.table[table_A]['BHID'].unique():
                af = groupA.get_group(i)['FROM'].values
                at = groupA.get_group(i)['TO'].values
                ai = groupA.get_group(i)['_ida'].values
            else: 
                af = []
                at = []
                ai = []
            if i in self.table[table_B]['BHID'].unique():
                bf = groupB.get_group(i)['FROM'].values
                bt = groupB.get_group(i)['TO'].values
                bi = groupB.get_group(i)['_idb'].values
            else: 
                bf = []
                bt = []
                bi = []


            # merge a on b
            zf,zt,zi = __intersect(af,at,ai, bf,bt)
            taf = taf + list(zf) # append to list
            tat = tat + list(zt)
            tai = tai + list(zi)
            ta_bhid = ta_bhid + [i for k in range(len(zf))]
             
            # merge b on a
            zf,zt,zi = __intersect(bf,bt,bi, af,at)
            tbf = tbf + list(zf) # append to list
            tbt = tbt + list(zt)
            tbi = tbi + list(zi)
            tb_bhid = tb_bhid + [i for k in range(len(zf))]

        #create dataframes of the two tables 
        b2a = pd.DataFrame({'BHID':ta_bhid, 'FROM':taf, 'TO':tat, '_ida':tai, })
        a2b = pd.DataFrame({'BHID':tb_bhid, 'FROM':tbf, 'TO':tbt, '_idb':tbi, })

        # the skeleton of the two tables with intervals intersected
        newtable = b2a.merge(a2b, on = ['BHID','FROM', 'TO'], how = 'outer').sort_values(by = ['BHID','FROM', 'TO'])

        # merge with existing data
        newtable=newtable.join(self.table[table_A], on='_ida', rsuffix='__tmp__')
        newtable=newtable.join(self.table[table_B], on='_idb', rsuffix='__tmp__')


        #clean if necessary
        if clean:
            newtable.drop(
               ['BHID__tmp__', 'FROM__tmp__','TO__tmp__','_ida__tmp__','_idb__tmp__',
                'xm','ym','zm','xb','yb','zb','xe','ye','ze',
                'xm__tmp__','ym__tmp__','zm__tmp__',
                'xb__tmp__','yb__tmp__','zb__tmp__',
                'xe__tmp__','ye__tmp__','ze__tmp__',
                ],
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
        """fix_zero_interval(str table_name, str new_table_name, bint overwrite = False, double tol = 0.01, bint clean = True, bint endhole=False, bint addgaps=True)

        Removes zero length intervals

        The function removes zero intervals and
        adds gaps if the gap length is longer than the tolerance.



        Parameters
        ----------
        table_name: name of the table

        new_table_name: name of the new table
            it may not exists at drillhole.table.keys()
        overwrite: default True.
            If new_table_name exists and overwrite == True the existing
            table will be overwritten.
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


        Example
        -------
        >>>
        >>> gap,overlap = mydrillhole.fix_zero_interval(table_name= 'assay',
                                                new_table_name = 'tmp',
                                                overwrite = False,
                                                tol = 0.01,
                                                clean = True,
                                                endhole = False,
                                                addgaps = True)
        >>>



        Note
        ----
        This function removes zero intervals and the calls the add_gaps
        function.


        TODO
        ----
         - [] Fix example section.

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



    cpdef collar2table(self, str table_name,
                       str new_table_name,
                       object collar_prop,
                       bint overwrite=False):
        """collar2table(str table_name, str new_table_name, list collar_prop, bint overwrite=False)

        Add collar properties to a table.


        Parameters
        ----------
        table_name : name of the table
        new_table_name : name of the output table
        collar_prop : list with property names in collar
            the property names may exist at drillhole.collar.columns
            if collar_prop= Null all the Collar properties will be copied
        overwrite : default False
            If new_table_name exists and overwrite == True the existing
            table will be overwrite.

        Example
        -------
        >>>
        >>> mydrillhole.collar2table(table_name = 'assay',
                                     new_table_name = 'assay',
                                     collar_prop = ['TYPE', 'COMMENT'],
                                     overwrite = True)
        >>>

        TODO
        ----
         - [] Fix example section.

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

        """downh_composite(str table_name, str variable_name, str new_table_name, double cint = 1, double minlen=-1, bint overwrite =False)

        Downhole composites one variable at the time


        Parameters
        ----------
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
            will be overwritten.

        Example
        -------
        >>>
        >>> mydrillhole.downh_composite(table_name = 'assay',
                                        variable_name = 'Au',
                                        new_table_name = 'cmp',
                                        cint = 1,
                                        minlen =-1,
                                        overwrite = False)
        >>>

        Note
        ---
        Undefined intervals in 'variable_name' will be excluded.

        To composite per lithology, you can filter out each lithology
        composite and then combine the composite tables.

        This algorithm starts from distance zero and no residual will
        be created at the top, for example, the if the first interval is
        [FROM=1.5, TO=2] it will be composited as [FROM=0, TO=2,_len=0.5]
        for 2 m composites. This produce consistent results if different
        tables are composited, making easier later table combination

        The accumulation output '_acum' can be used to produce accumulated
        variables, to identify number of intervals in each composite or
        to track categorical variables (coded as number) compositing.

        TODO
        ----
         - [] Fix example section.
         - [] Add multiple variables
         - [] Add domain code

        """

        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)

        #check that the variable is in the table and the variable type
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

            nf, nt, nlen, nvar, nacum= __composite1dh(ifrom= group.get_group(i)['FROM'].values,
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

    cpdef key_composite(self,  str table_name,
                          str key_name,
                          str variable_name,
                          str new_table_name,
                          double tol=0.1,
                          bint overwrite =False):

        """downh_composite(str table_name, str key_name, str variable_name, str new_table_name, double tol=0.001, bint overwrite =False)

        Downhole composites one variable at the time


        Parameters
        ----------
        table_name: name of the table with drillhole intervals
            it must be an existing table at drillhole.table.keys()
        key_name: name of the key in "table_name" table
            it must be an existing field at drillhole.table[table_name].columns
        variable_name: name of the variable in "table_name" table
            it must be an existing field at drillhole.table[table_name].columns
        new_table_name: name of the new table with composites
            it can be a new name, no in drillhole.table.keys()
            or a table name on drillhole.table.keys() if overwrite =True
        tol: default 0.1, tolerance to ignore gaps
        overwrite: default True.
            If the table exist and overwrite = True the existing table
            will be overwritten.

        Example
        -------
        >>>



        TODO
        ----
         - [] Fix example section.
         - [] Add accumulator and count variable

        """



        cdef double l, vl, f, t
        cdef int i, n

        # check that the table is not in the database
        if overwrite==False:
            assert new_table_name not in self.table.keys(), 'The table {} already exist, use overwrite = True to rewrite'.format(new_table_name)

        #check that the variable is in the table and the variable type
        assert variable_name in self.table[table_name].columns, 'Error: The variable {} do not exists in the input table'.format(variable_name)
        assert np.dtype(self.table[table_name][variable_name])==np.float64, 'Error: The variable {} is not type float64'.format(variable_name)
        assert key_name in self.table[table_name].columns, 'Error: The key variable {} do not exists in the input table'.format(key_name)

        # may not work with nrows < 2
        assert self.table[table_name].shape[0]>2, 'Error: This function requires tables with more than one row'

        # sort table
        self.table[table_name].sort_values(by=['BHID', 'FROM'], inplace=True)


        #create input arrays
        ibhid =  self.table[table_name]['BHID'].values
        ifrom = self.table[table_name]['FROM'].values
        ito =   self.table[table_name]['TO'].values
        ikey =  self.table[table_name][key_name].values
        ivar =  self.table[table_name][variable_name].values

        #create output arrays
        nrows = self.table[table_name].shape[0]
        obhid = np.empty(nrows, dtype = self.table[table_name]['BHID'].dtype )
        ofrom = np.empty(nrows, dtype = self.table[table_name]['FROM'].dtype)
        oto = np.empty(nrows, dtype = self.table[table_name]['TO'].dtype)
        okey = np.empty(nrows, dtype = self.table[table_name][key_name].dtype)
        ovar = np.empty(nrows, dtype = self.table[table_name][variable_name].dtype)


        #composite by key
        vl = ivar[0]*(ito[0]-ifrom[0])
        l = ito[0]-ifrom[0]
        f = ifrom[0]
        t = ito[0]
        key = ikey[0]
        bhid = ibhid[0]
        n = 0

        for i in range(1,nrows):

            if ibhid[i]==ibhid[i-1] and ikey[i]==ikey[i-1] and ifrom[i]-ito[i-1]<=tol:
                vl += ivar[i]*(ito[i]-ifrom[i])
                l  +=  ito[i]-ifrom[i]
                t   =  ito[i]

            else:

                obhid[n] = bhid
                ofrom[n] = f
                oto[n] = t
                okey[n] = key
                ovar[n] = vl/l

                # Update next iteration
                n+=1
                vl = ivar[i]*(ito[i]-ifrom[i])
                l = ito[i]-ifrom[i]
                f = ifrom[i]
                t = ito[i]
                key = ikey[i]
                bhid = ibhid[i]

        #take right action with las interval in the array

        # The nth interval ==  nth-1 interval, add interval to composit and close
        if ibhid[nrows-1]==ibhid[nrows-2] and ikey[nrows-1]==ikey[nrows-2] and ifrom[nrows-1]-ito[nrows-2]<=tol:

            obhid[n] = bhid
            ofrom[n] = f
            oto[n] = t
            okey[n] = key
            ovar[n] = vl/l
            n += 1


        # The nth interval !=  nth-1 interval so, add this single interval and close
        else:
            # Update next iteration
            obhid[n] = bhid
            ofrom[n] = f
            oto[n] = t
            okey[n] = key
            ovar[n] = vl/l
            n += 1

        #create new table with gaps (only with fields )
        newtable=pd.DataFrame({'BHID':obhid[:n], 'FROM':ofrom[:n],'TO':oto[:n], key_name:okey[:n], variable_name:ovar[:n]})

        #add table to the class
        self.addtable(newtable,new_table_name,overwrite)



    cpdef bench_composite(self, str table_name,
                         double zmax,
                         double bench,
                         double tol=0.01):
        """bench_composite(str table_name, double zmax, double bench, double tol=0.01)

        This function is not implemented yet.
        """


        print 'we are working on that'


    #-------------------------------------------------------------------
    #       VTK export
    #-------------------------------------------------------------------
    cpdef intervals2vtk(self, str table_name, str filename):
        """ intervals2vtk(str table_name, str filename)

        Export desurveyed drillhole table to vtk lines. Endpoints
        are required.


        Parameters
        ----------
        table_name : str
        filename : str
            This is the absolute or relative file path and name.

        See Also
        --------
        pygslib.vtktool


        Note
        ----
        To export to VTK points see pygslib.vtktool

        Examples
        --------
        >>>
        >>> mydrillhole.export_core_vtk_line(table_name = 'assay',
                                            filename = 'assay_line.vtk')
        >>>

        TODO
        ----
         - [] Fix example section

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

        cdef int l, n, dlen
        cdef np.ndarray[double, ndim=1] xb
        cdef np.ndarray[double, ndim=1] yb
        cdef np.ndarray[double, ndim=1] zb
        cdef np.ndarray[double, ndim=1] xe
        cdef np.ndarray[double, ndim=1] ye
        cdef np.ndarray[double, ndim=1] ze


        try:

            # this will work if all coordinates defined

            assert  np.isfinite(self.table[table_name]['xb'].values).all(), "no finite coordinates at xb, please filter out nonfinite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['yb'].values).all(), "no finite coordinates at yb, please filter out nonfinite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['zb'].values).all(), "no finite coordinates at zb, please filter out nonfinite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['xe'].values).all(), "no finite coordinates at xe, please filter out nonfinite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['ye'].values).all(), "no finite coordinates at ye, please filter out nonfinite coordinates and try again"
            assert  np.isfinite(self.table[table_name]['ze'].values).all(), "no finite coordinates at ze, please filter out nonfinite coordinates and try again"

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

            # if ther are undefined values at yb,zb, xe,ye,ze the output will be wrong, here we check this
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

            # copy data
            if  dtype==np.int8 or dtype==np.int16 or dtype==np.int32 or dtype==np.int64 or dtype==np.float16 or dtype==np.float32 or dtype==np.float64:
                vtkfields[i]= vtk.util.numpy_support.numpy_to_vtk(self.table[table_name][i].values)
                vtkfields[i].SetName(i)
                vtkfields[i].SetNumberOfComponents(1)
            else:
                # this is fos string array. Not optimized...
                vtkfields[i]= vtk.vtkStringArray()
                vtkfields[i].SetName(i)
                vtkfields[i].SetNumberOfComponents(1)
                vtkfields[i].SetNumberOfTuples(dlen)
                for l in range(dlen):
                    vtkfields[i].SetValue(l,str(self.table[table_name][i][l]))

        # now we create a set of vtk points
        points= vtk.vtkPoints()
        npoints = dlen*2
        points.SetNumberOfPoints(npoints)

        # now we create a set of lines representing the cores and
        # a line container (a cell array)
        line = vtk.vtkLine()
        lines = vtk.vtkCellArray()


        # TODO: Optimize this for... replace points.InsertNextPoint with C array.


        # populate this data
        n=-1
        for l in range(dlen):

            n=n+1
            points.SetPoint(n,xb[l], yb[l], zb[l])
            line.GetPointIds().SetId(0,n)
            n=n+1
            points.SetPoint(n,xe[l], ye[l], ze[l])
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
        pygslib.vtktools.SavePolydata(linesPolyData, filename)
