# simple sum function
import numpy as np
cpdef zadd(double[:] a,double[:] b):
      cdef size_t i;   
      c = np.zeros((a.shape[0]), dtype = np.float64)
      cdef double[:] c_view = c    
      for i in range(a.shape[0]):
            c[i] = a[i] + b[i]
      return c
