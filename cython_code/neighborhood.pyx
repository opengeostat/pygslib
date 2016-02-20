'''
PyGSLIB Neighborhood, Module to handle search on data.  

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

from scipy import spatial
import pygslib as gslib    
import numpy as np
cimport numpy as np
import pandas as pd

cdef class Neighborhood:
    
    cdef readonly:
        # data and target array
        # pass data by reference
        np.ndarray target_x
        np.ndarray target_y 
        np.ndarray target_z 
        np.ndarray data_x 
        np.ndarray data_y 
        np.ndarray data_z 
        
        np.ndarray target_x_r
        np.ndarray target_y_r
        np.ndarray target_z_r 
        np.ndarray data_x_r
        np.ndarray data_y_r
        np.ndarray data_z_r
        
        object target_properies
        object data_properies 
    
        
        # search parameters
        double axis_x 
        double axis_y 
        double axis_z 
        
        double ang1
        double ang2
        double ang3
        
        double anis1 
        double anis2 
        
        #get some data properties
        int nt 
        int nd
        int npt
        int npd
        
        #ids
        np.ndarray target_ids
        np.ndarray data_ids 
        
        int current_target
        
        np.ndarray selected_data 
        np.ndarray distance_data2target 
        
        
        # parameters for rotation
        object parameters_data      
        object parameters_target          
    
        # the tree structure for search
        object tree 
    
    
    def __cinit__(self, 
                 np.ndarray[double, ndim=1] target_x,
                 np.ndarray[double, ndim=1] target_y, 
                 np.ndarray[double, ndim=1] target_z,
                 np.ndarray[double, ndim=1] data_x,
                 np.ndarray[double, ndim=1] data_y,
                 np.ndarray[double, ndim=1] data_z,
                 object target_properies,
                 object data_properies, 
                 object angles,
                 object axis,
                 long leafsize=1000000000):
        
        # some sheck
        assert target_z.shape[0]==target_y.shape[0]==target_x.shape[0], 'Error: target coordinates with different num of rows'
        assert data_z.shape[0]==data_y.shape[0]==data_x.shape[0], 'Error: data coordinates with different num of rows'
        assert target_properies.shape[0]==target_x.shape[0], 'Error: target coordinates and properies arrays with different num of rows'
        assert data_properies.shape[0]==data_x.shape[0], 'Error: data coordinates and properies arrays with different num of rows'
        assert axis[0] > 0, 'search axis[0]== 0, use all axis > 0'
        assert axis[1] > 0, 'search axis[1]== 0, use all axis > 0'
        assert axis[2] > 0, 'search axis[2]== 0, use all axis > 0'
        
        
        
        # data and target array
        # pass data by reference
        self.target_x = target_x
        self.target_y = target_y
        self.target_z = target_z
        self.data_x = data_x
        self.data_y = data_y
        self.data_z = data_z
        self.target_properies = target_properies
        self.data_properies = data_properies
    
        
        # search parameters
        self.axis_x = axis[0]
        self.axis_y = axis[1]
        self.axis_z = axis[1]
        
        self.ang1 = angles[0]
        self.ang2 = angles[1]
        self.ang3 = angles[1]
        
        self.anis1 = self.axis_y / self.axis_x
        self.anis2 = self.axis_z / self.axis_x
        

        #get some data properties
        self.nt = self.target_x.shape[0] # num of targets
        self.nd = self.data_x.shape[0]
        
        if self.target_properies.ndim > 1:
            self.npt = self.target_properies.shape[1] # num of targets
        else: 
            self.npt = 1
            
        if self.data_properies.ndim > 1:
            self.npd = self.data_properies.shape[1]
        else:
            self.npd = 1
        
        self.target_ids = np.arange(self.nt, dtype=int)
        self.data_ids = np.arange(self.nd, dtype=int)
        
        self.current_target = 0
        
        self.selected_data = np.zeros(self.nd, dtype=int)
        self.distance_data2target = np.zeros(self.nd, dtype=float)
        
        
        # prepare a new dataset with transformed coordinates to an isotropic system
        # first we define the parameters for rotation
        self.parameters_data = { 
                            'x'      :  self.data_x,
                            'y'      :  self.data_y,
                            'z'      :  self.data_z,
                            'x0'     :  0,           # new X origin of coordinate , 'f' 
                            'y0'     :  0,           # new Y origin of coordinate , 'f'
                            'z0'     :  0,           # new Z origin of coordinate , 'f'
                            'ang1'   :  self.ang1,   # Z  Rotation angle, 'f' 
                            'ang2'   :  self.ang2,   # X  Rotation angle, 'f'
                            'ang3'   :  self.ang3,   # Y  Rotation angle, 'f'
                            'anis1'  :  self.anis1,  # Y cell anisotropy, 'f' 
                            'anis2'  :  self.anis1,  # Z cell anisotropy, 'f' 
                            'invert' :  0}           # 0 do rotation, <> 0 invert rotation, 'i'

        self.parameters_target = { 
                            'x'      :  self.target_x,
                            'y'      :  self.target_y,
                            'z'      :  self.target_z,
                            'x0'     :  0,           # new X origin of coordinate , 'f' 
                            'y0'     :  0,           # new Y origin of coordinate , 'f'
                            'z0'     :  0,           # new Z origin of coordinate , 'f'
                            'ang1'   :  self.ang1,   # Z  Rotation angle, 'f' 
                            'ang2'   :  self.ang2,   # X  Rotation angle, 'f'
                            'ang3'   :  self.ang3,   # Y  Rotation angle, 'f'
                            'anis1'  :  self.anis1,  # Y cell anisotropy, 'f' 
                            'anis2'  :  self.anis1,  # Z cell anisotropy, 'f' 
                            'invert' :  0}           # 0 do rotation, <> 0 invert rotation, 'i'
        
        # then we rotate
        self.data_x_r,self.data_y_r,self.data_z_r = gslib.rotscale(self.parameters_data)
        self.target_x_r,self.target_y_r,self.target_z_r = gslib.rotscale(self.parameters_target)
        
        
        # we create a tree structure for search
        # warning: This array is not copied, and so modifying this data will result in bogus results.
        
        # important leafsize realy large is required for some data sets: see 
        # http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.cKDTree.html
        
        self.tree = spatial.cKDTree(pd.DataFrame({'x':self.data_x_r,
                                                 'y':self.data_y_r,
                                                 'z':self.data_z_r}), leafsize=leafsize)
        
        
    def __iter__(self):
        return self

    cpdef next(self):
        if self.current_target < self.nt:
            self.current_target += 1
            
            # return ???
            return 0  # TODO
        else:
            raise StopIteration()
    
    cpdef get_data(self, int target=0, double scale=1., int k=1): 
        
        
        # check bounds
        assert target<self.nt, 'target = {} out of bound, use target < {}'.format(target, self.nt)

        #set current tarhet for external reference
        self.current_target=target 
        
        #select grup of points with the tree
        self.distance_data2target, self.selected_data = self.tree.query([self.target_x_r[target],
                              self.target_y_r[target],
                              self.target_z_r[target]], k=k, distance_upper_bound=scale*self.axis_x)
        
        oks = np.isfinite(self.distance_data2target)
        self.distance_data2target = self.distance_data2target[oks]
        self.selected_data = self.selected_data[oks]
                              
    
    cpdef apply_interpolator(self, object interpolator, object parameters={}, 
                             object dtypes=[float], double scale=1., object k=1):
        
        
        cdef:
            int i
            int noutputs
            int nt
        
        #define outputs
        noutputs=len(dtypes)
        output=[]
        for i in range(noutputs): 
            output.append(np.zeros(self.nt, dtype=dtypes[i]))

        
        # apply function
        nt = self.nt
        for i in range(self.nt):
            
                
            # select data in neighborhood i (selection is applied internally on self.selected_data and self.distance_data2target)
            self.get_data(target=i, scale=scale, k=k)

            if len(self.distance_data2target)>0:

                # apply funtion 
                # warning, this function may defined to access data from self.selected_data,  self.distance_data2target, etc.
                tmpout = interpolator(self, **parameters)

                # populate results
                for j in range(noutputs):
                    output[j][i]= tmpout[j]
        
        #return results
        return output


    cpdef test(self):
        
        
        cdef:
            int i
            int nt
        
        #define outputs
        output=np.empty(self.nt)
    
        # apply function
        nt = self.nt
        for i in range(self.nt):
            
            # select data in neighborhood i (selection is applied internally on self.selected_data and self.distance_data2target)
            #self.get_data(target=i, scale=1., k=10)
            
            if len(self.distance_data2target)>0:
                pass
        
        #return results
        return output
