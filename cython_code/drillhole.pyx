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
#  General functions
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
    cdef float x
    cdef float y
    cdef float z

    # internal
    cdef float razm
    cdef float rdip
    cdef float DEG2RAD

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
    cdef float azm
    cdef float dip


    # internal
    cdef float razm
    cdef float rdip
    cdef float RAD2DEG
    cdef float pi

    RAD2DEG=180.0/3.141592654
    pi = 3.141592654

    if x!=0. and y!= 0.: 
        azm= atan2(x,y)
        if azm<0:  
            azm= azm + pi*2
        azm = azm * RAD2DEG
    else: 
        azm = 0


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
    cdef float azm
    cdef float dip

    # internal
    cdef float x1
    cdef float y1
    cdef float z1
    cdef float x2
    cdef float y2
    cdef float z2
    cdef float x
    cdef float y
    cdef float z


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
    cdef float dz
    cdef float dn
    cdef float de 


    # internal 
    cdef float i1
    cdef float a1
    cdef float i2
    cdef float a2
    cdef float DEG2RAD
    cdef float rf
    cdef float dl

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
               float lpt):

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


    # output (anges at begin, mid and end interval)
    cdef float azt
    cdef float dipt
    cdef float xt
    cdef float yt
    cdef float zt

    # internal
    cdef int i
    cdef int j
    cdef float a,
    cdef float b
    cdef float azm1
    cdef float dip1
    cdef float azm2
    cdef float dip2
    cdef float len12
    cdef float d1
    cdef float EPSLON=1.0e-4
    cdef float xa
    cdef float ya
    cdef float za
    cdef float xb
    cdef float yb
    cdef float zb
    cdef float dz
    cdef float dn
    cdef float de

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
            # desurvey interval at zurvey... table

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
        xt = xa+de
        yt = ya+dn
        zt = za-dz 
        
        warnings.warn('\n point beyond the last survey point at %s' % indes)
    else:
        warnings.warn('\n not interval found at survey, at %s' % indes)

    return   azt, dipt, xt, yt, zt

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
        
        #check we have the rigth naming in survey columns 
        assert 'BHID' in survey.columns, "survey don't have BHID column"
        assert 'AT' in survey.columns, "survey don't have AT column"
        assert 'AZ' in survey.columns, "survey don't have AZ column"
        assert 'DIP' in survey.columns, "survey don't have DIP column"
        
        
        self.collar = collar
        self.survey = survey
        self.table = {}
        
           
    
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
            self.table[table_name]=table
        else:
            if overwrite == True: 
                self.table[table_name]=table
            else:
                raise NameError('Table %s already exist, use overwrite = True to overwrite' % table_name)
            
        
    cpdef validate(self):
        """
        validate()
        
        Runs a set of basic validations on survey and collar 
        consisting in: 
        - check existence of Null values
        - check survey without values AT=0
        - check that coordinates and direction angles dtypes are float64
        
        Notes
        -----
        a) You may run this validation before doing desurvey
        b) Ony few minimum validation are implemented for now
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
        
        #check collar
        # null values in collar
        if self.collar['BHID'].hasnans():
            raise NameError('Non defined BHID in collar table')
        if self.collar['XCOLLAR'].hasnans():
            raise NameError('Non defined XCOLLAR in collar table')
        if self.collar['YCOLLAR'].hasnans():
            raise NameError('Non defined YCOLLAR in collar table')
        if self.collar['ZCOLLAR'].hasnans():
            raise NameError('Non defined ZCOLLAR in collar table')
        if self.collar['XCOLLAR'].dtypes!='float64':
            raise NameError('XCOLLAR in collar table != float64')            
        if self.collar['YCOLLAR'].dtypes!='float64':
            raise NameError('YCOLLAR in collar table != float64')
        if self.collar['ZCOLLAR'].dtypes!='float64':
            raise NameError('ZCOLLAR in collar table != float64')

        #check SURVEY
        # null values in survey
        if self.survey['BHID'].hasnans():
            raise NameError('Non defined BHID in survey table')
        if self.survey['AT'].hasnans():
            raise NameError('Non defined AT in survey table')
        if self.survey['AZ'].hasnans():
            raise NameError('Non defined AZ in survey table')
        if self.survey['DIP'].hasnans():
            raise NameError('Non defined DIP in survey table')
        if self.survey['AT'].dtypes!='float64':
            raise NameError('AT in survey table != float64')            
        if self.survey['DIP'].dtypes!='float64':
            raise NameError('DIP in survey table != float64')
        if self.survey['AZ'].dtypes!='float64':
            raise NameError('AZ in survey table != float64')
        
        #check survey without values: at=0
        self.survey.sort(columns=['BHID','AT'], inplace=True)
        error = self.__checkAt0(self.survey['BHID'].values, self.survey['AT'].values)
        if error>-1:
            raise NameError('Firts inteval AT!=0 at survey table, positiom %d' %error) 
            
        return None
        # TODO: check table relationship
        
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
        # null values in table bhid
        if self.table[table_name]['BHID'].hasnans():
            raise NameError('Non defined BHID in %s' % table_name)
        # null values in From/To
        if self.table[table_name]['FROM'].hasnans():
            raise NameError('Non defined FROM in %s' % table_name)
        if self.table[table_name]['TO'].hasnans():
            raise NameError('Non defined TO in %s' % table_name)
        if self.table[table_name]['FROM'].dtypes!='float64':
            raise NameError('FROM in table %s != float64' % table_name)
        if self.table[table_name]['TO'].dtypes!='float64':
            raise NameError('TO in table %s != float64' % table_name)

        # TODO: check overlaps and table relationship

    
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



    cpdef desurvey(self, str table_name, bint endpoints=False):
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
        
        
        cdef np.ndarray[object, ndim=1] idc = self.collar['BHID'].values
        cdef np.ndarray[double, ndim=1] xc = self.collar['XCOLLAR'].values
        cdef np.ndarray[double, ndim=1] yc = self.collar['YCOLLAR'].values
        cdef np.ndarray[double, ndim=1] zc = self.collar['ZCOLLAR'].values
        cdef np.ndarray[object, ndim=1] ids = self.survey['BHID'].values
        cdef np.ndarray[double, ndim=1] ats = self.survey['AT'].values
        cdef np.ndarray[double, ndim=1] azs = self.survey['AZ'].values
        cdef np.ndarray[double, ndim=1] dips = self.survey['DIP'].values
        cdef np.ndarray[object, ndim=1] idt =self.table[table_name]['BHID'].values
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
              
        indbt = 0
        indet = 0 
        inds = 0
        indt = 0   
        for jc in range(nc):
            indbs = -1
            indes = -1
            # get first index of the collar jc in survey 
            for js in range(inds, ns):
                # find index beging and end for the actual collar
                if idc[jc]==ids[js]:
                    inds = js
                    indbs = js
                    break
            
            # get las index of the collar jc in survey 
            for js in range(inds, ns):
                # find index beging and end for the actual collar
                if idc[jc]!=ids[js]:
                    break
                else: 
                    inds = js
                    indes = js

            if indbs==-1 or indes==-1:
                # do not desurvey this drillhole
                warnings.warn('! collar without survey, table not desurveyed')
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
                    desurv1dh(indbs,indes,ats,azs,dips,xc[jc],yc[jc],zc[jc],mid)
                    
                    if endpoints==True:
                        tmpaz,tmpdip,xbt[jt],ybt[jt],zbt[jt] = \
                        desurv1dh(indbs,indes,ats,azs,dips,xc[jc],yc[jc],zc[jc],fromt[jt])
                        tmpaz,tmpdip,xet[jt],yet[jt],zet[jt] = \
                        desurv1dh(indbs,indes,ats,azs,dips,xc[jc],yc[jc],zc[jc],tot[jt])  
        
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
            
