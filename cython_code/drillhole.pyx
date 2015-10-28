
import pygslib as gslib
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
%load_ext Cython
%%cython
import numpy as np
cimport numpy as np
def txt2intID(np.ndarray[str, ndim=1] cTextID, np.ndarray[str, ndim=1] tTextID):
    """
    arange(np.ndarray[str, ndim=1] cTextID, np.ndarray[str, ndim=1] tTextID)
    
    Return two numpy 1D arrays representing the integer ID in a collar 
    and any other acompaning table. 
    
    Imput values may be sorted with the same criteria. ID in assay 
    but not in collar will produce undefined values. 
    
    
    Parameters
    ----------
    cTextID : 1D numpy array with type str 
        This is the collar BHID and may be unique and sorted. 
    tTextID : 1D numpy array with type str 
        This is the Table  BHID and may be unique and sorted. 
    
    Returns
    -------
    cBHID : ndarray (1D, dtype=int) 
        Array of BHID at collar table.
    
    tBHID : ndarray (1D, dtype=int) 
        Array of BHID at table matching collar BHIDs.
    
    See Also
    --------
    
    Examples
    --------
    >>> # Adding integer BHID to collar and survey table in pandas DataFrames
    >>> sBHID,cBHID = txt2intID(collar['HOLE-ID'].values, survey['HOLE-ID'].values)
    >>> collar['BHID']= cBHID
    >>> survey['BHID']= sBHID
    >>> # Adding integer BHS to assay Table
    >>> aBHID,cBHID = txt2intID(collar['HOLE-ID'].values, assays['HOLE-ID'].values)
    >>> assays['BHID']= aBHID
    >>> # check results, example NaN values (assay without collar)
    >>> assays[np.isnan(assays['BHID'])]
    """
    
    # the data is assumed sorted
    
    cdef int sloc
    cdef int i
    cdef int j
    cdef int nc = cTextID.shape[0]
    cdef int nt = tTextID.shape[0]
    cdef np.ndarray cBHID = np.zeros([nc], dtype=int)
    cdef np.ndarray tBHID = np.empty([nt], dtype=int)
    
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
    return tBHID,cBHID
