# -*- coding: utf-8 -*-
#!/usr/bin/env python

# using naming on http://www.gslib.com/gslib_help/programs.html
import subprocess
import copy
import pandas as pd
import pygslib
import numpy as np
import os

__varmap_par = \
"""                  Parameters for VARMAP
                  *********************

START OF PARAMETERS:
{datafl}                          -file with data
{nvar} {ivar_}               - number of variables: column numbers
{tmin} {tmax}               - trimming limits
{igrid}                     - 1=regular grid, 0=scattered values
{nx} {ny} {nz}              - if igrid=1: nx, ny, nz
{xsiz} {ysiz} {zsiz}        -  xsiz, ysiz, zsiz
{icolx} {icoly} {icolz}     - if igrid=0: columns for x, y, z coordinates
{outfl}                     -file for variogram grid output
{nxlag} {nylag} {nzlag}     -nxlag, nylag, nzlag
{dxlag} {dylag} {dzlag}     -dxlag, dylag, dzlag
{minpairs}                  -minimum number of pairs
{standardize}               -standardize sill? (0=no, 1=yes)
{nvarg}                     -number of variograms
{ivpar_}                    -tail, head, variogram type (array with shape [nvarg,4])


 vg type 1 = traditional semivariogram
         2 = traditional cross semivariogram
         3 = covariance
         4 = correlogram
         5 = general relative semivariogram
         6 = pairwise relative semivariogram
         7 = semivariogram of logarithms
         8 = semimadogram
         9 = indicator semivariogram - continuous
         10= indicator semivariogram - categorical

cut[i] is only required if ivtype[i] == 9 or == 10

"""

def varmap(parameters, gslib_path = None, silent = False, xorg=0., yorg=0., zorg=0.):
    """varmap(parameters, gslib_path = None)

    Funtion to calculate variogram maps (grid) using varmap.exe external
    gslib program.

    Parameters
    ----------
    parameters : dict
        dictionary with parameters
    gslib_path : string (default None)
        absolute or relative path to gslib excecutable programs
    silent: boolean
        if true external GSLIB stdout text is printed
    xorg, yorg, zorg: floats (default 0.)
        origin of coordinated of the variogram map in the vtkImageData output
    Returns
    ------
    (pandas.DataFrame, vtkImageData) with variogram map results

    Example
    --------
    TODO:

    Notes
    ------
    The dictionary with parameters may be as follows::


        parameters = {
            'datafl' : str or, None, or numpy,     # path to file, or none (to use '_xxx_.in') or numpy array (with columns [x,y])
            'ivar'   : 1D array of int,            # variables column numbers to be used in ivtail and ivhead,
            'tmin'   : float,                      # trimming limits min and max (raws out of this range will be ignored)
            'tmax'   : float,
            'igrid'   : int,                        # 1=regular grid, 0=scattered values
            'nx':int,'ny':int,'nz':int,             # if igrid=1: nx, ny, nz
            'xsiz':float,'ysiz':float,'zsiz':float, # if igrid=1: xsiz, ysiz, zsiz
            'icolx':int,'icoly':int,'icolz':int,    # if igrid=0: columns for x, y, z coordinates
            'outfl': str or None,                   # path to the output file or None (to use '_xxx_.out')
            'nxlag':int,'nylag':int,'nzlag':int,    # nxlag, nylag, nzlag
            'dxlag':float,'dylag':float,'dzlag':float,    # dxlag, dylag, dzlag
            'minpairs': int,                        # minimum number of pairs
            'standardize': int,                     # standardize sill? (0=no, 1=yes)
            'ivpar': 2D array of int                # tail, head, variogram type, and cut (with shape [nvarg,4])
            }


            ivtype   1 = traditional semivariogram
                     2 = traditional cross semivariogram
                     3 = covariance
                     4 = correlogram
                     5 = general relative semivariogram
                     6 = pairwise relative semivariogram
                     7 = semivariogram of logarithms
                     8 = semimadogram
                     9 = indicator semivariogram - continuous
                     10= indicator semivariogram - categorical

    see http://www.gslib.com/gslib_help/varmap.html for more information

    """

    if gslib_path is None:
        if os.name == "posix":
            gslib_path = '~/gslib/varmap'
        else:
            gslib_path = 'c:\\gslib\\varmap.exe'

    mypar = copy.deepcopy(parameters)

    # handle the case where input is an array an not a file
    if isinstance(parameters['datafl'], np.ndarray):
        mypar['datafl']='_xxx_.in'
        if parameters['datafl'].ndim<2:
            parameters['datafl']= parameters['datafl'].reshape([parameters['datafl'].shape[0],1])
        if mypar['igrid']== 1:
            mypar['ivar'] = np.arange(parameters['datafl'].shape[1])+1
        else:
            mypar['ivar'] = np.arange(parameters['datafl'].shape[1]-3)+1
            mypar['icolx']=1
            mypar['icoly']=2
            mypar['icolz']=3
        with open('_xxx_.in',"w") as f:
            f.write('temp file '+'\n')
            f.write('{}'.format(parameters['datafl'].shape[1])+'\n')
            if mypar['igrid']==1:
                for i in range(parameters['datafl'].shape[1]):
                    	f.write('v{}\n'.format(i+1))
            else:
                f.write('x\ny\nz\n')
                for i in range(3,parameters['datafl'].shape[1]):
                    	f.write('v{}\n'.format(i-2))
            np.savetxt(f,parameters['datafl'])

    elif parameters['datafl'] is None:
        mypar['datafl']='_xxx_.in'


    # some updates to ensure the parameterfile is good
    if mypar['igrid']==0:
        mypar['nx']= 0
        mypar['ny']= 0
        mypar['nz']= 0
        mypar['xsiz']= 0
        mypar['ysiz']= 0
        mypar['zsiz']= 0

    if mypar['igrid']==1:
        mypar['icolx']= 0
        mypar['icoly']= 0
        mypar['icolz']= 0

    if mypar['outfl'] is None:
        mypar['outfl'] = '_xxx_.out'

    ivpar = np.array (mypar['ivpar'])
    assert (ivpar.shape[1]==4)
    assert (set(ivpar[:,0]).issubset(set(mypar['ivar'])))  # head variable
    assert (set(ivpar[:,1]).issubset(set(mypar['ivar'])))  # tail variable
    assert (set(ivpar[:,2]).issubset(set([1,2,3,4,5,6,7,8,9,10]))) # ivtype
    for i in range(ivpar.shape[0]):
        if ivpar[i,2]<9:
            ivpar[i,3] = None
        else:
            if ivpar[i,3]==None:
                raise NameError('gslib varmap Error inparameter file: cut[{}]=None'.format(i))


    # prepare parameter file and save it
    mypar['nvar'] = len(mypar['ivar'])
    mypar['ivar_'] = ' '.join(map(str, mypar['ivar'])) # array to string
    mypar['nvarg'] = ivpar.shape[0]
    mypar['ivpar_'] = pd.DataFrame.to_string(pd.DataFrame(ivpar),index= False, header=False) # array to string

    par = __varmap_par.format(**mypar)
    print (par)
    fpar ='_xxx_.par'
    with open(fpar,"w") as f:
        f.write(par)

    # call pygslib
    # this construction can be used in a loop for parallel execution
    p=subprocess.Popen([gslib_path, fpar],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    result = p.returncode
    p.wait()



    if p.returncode!=0:
        raise NameError('gslib declus NameError' + str(stderr.decode('utf-8')))

    if ~silent:
        try:
            print (stdout.decode('utf-8'))
        except:
            print (stdout)

    # return results as panndas array
    nxg=mypar['nxlag']*2+1
    nyg=mypar['nylag']*2+1
    nzg=mypar['nzlag']*2+1
    result = pygslib.gslib.read_gslib_file(mypar['outfl'])
    result['ID'] = np.repeat(np.arange(ivpar.shape[0]), nxg*nyg*nzg) # add variogram id
    # prepare array of vtk obfect with variogram map
    vmaps = []
    for i in range(mypar['nvarg']):
        vmaps.append(pygslib.vtktools.grid2vtkImageData(nx=nxg, ny=nyg, nz=nzg,
                                            xorg=xorg, yorg=yorg, zorg=zorg,
                                            dx=mypar['dxlag'],
                                            dy=mypar['dylag'],
                                            dz=mypar['dzlag'],
                                            cell_data=result[result['ID']==i].to_dict(orient='list')))



    return result, vmaps
