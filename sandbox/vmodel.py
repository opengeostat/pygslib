# -*- coding: utf-8 -*-
#!/usr/bin/env python

# using naming on http://www.gslib.com/gslib_help/programs.html
import subprocess
import copy
import pandas as pd
import pygslib
import numpy as np
import os
import matplotlib.pyplot as plt

__vmodel_par = \
"""                  Parameters for VMODEL
                  *********************

START OF PARAMETERS:
{outfl}                           - file for variogram output
{ndir} {nlag}                     - number of directions and lags
{ivdir_}                          - azm, dip, lag distance (array with shape [ndir,3])
{nst} {c0}                        - nst, nugget effect
{vst_}
"""

def vmodel(parameters, gslib_path = None, silent = False):
    """vmodel(parameters, gslib_path = None)

    Funtion to calculate variogram models using
    the vmodel.exe external gslib program.

    Parameters
    ----------
    parameters : dict
        dictionary with parameters
    gslib_path : string (default None)
        absolute or relative path to gslib excecutable programs
    silent: boolean
        if true external GSLIB stdout text is printed

    Returns
    ------
    pandas.DataFrame with variograms

    Example
    --------
    TODO:

    Notes
    ------
    The dictionary with parameters may be as follows::


        parameters = {
            'outfl'  : str or None,                  # path to the output file or None (to use '_xxx_.out')
            'nlag'   : int,                          # number of lags
            'ivdir'  : 2D array of floats,           # azm,dip,lag distance (array with shape [ndir,3])
            'c0'     : float,                        # nugget effect
            'it'     : int,                          # list of structure types (array with shape [nst])
            'cc'     : float,                        # list of structure variances (array with shape [nst])
            'ang1'   : float,                        # list of structure 1st rotation (array with shape [nst])
            'ang2'   : float,                        # list of structure 2nd rotation (array with shape [nst])
            'ang3'   : float,                        # list of structure 3rd rotation (array with shape [nst])
            'a_hmax' : float,                        # list of structure maximum horizontal ranges (array with shape [nst])
            'a_hmin' : float,                        # list of structure minimum horizontal ranges (array with shape [nst])
            'a_vert' : float}                        # list of structure vertical ranges (array with shape [nst])

    the parameters nst (number of structures) and ndir (number of directions)
    are theduced from arrays ivdir and it

    see http://www.gslib.com/gslib_help/vmodel.html for more information

    """

    if gslib_path is None:
        if os.name == "posix":
            gslib_path = '~/gslib/vmodel'
        else:
            gslib_path = 'c:\\gslib\\vmodel.exe'

    mypar = copy.deepcopy(parameters)

    if mypar['outfl'] is None:
        mypar['outfl'] = '_xxx_.out'

    # handle parameter arrays
    ivdir = np.array (mypar['ivdir'])
    it = np.array (mypar['it'], dtype = int)
    cc = np.array (mypar['cc'])
    ang1 = np.array (mypar['ang1'])
    ang2 = np.array (mypar['ang2'])
    ang3 = np.array (mypar['ang3'])
    a_hmax = np.array (mypar['a_hmax'])
    a_hmin = np.array (mypar['a_hmin'])
    a_vert = np.array (mypar['a_vert'])

    assert (ivdir.ndim==2)
    assert (it.ndim==cc.ndim==ang1.ndim==ang2.ndim==ang3.ndim==a_hmax.ndim==a_hmin.ndim==a_vert.ndim==1)
    assert (it.shape[0]==cc.shape[0]==ang1.shape[0]==ang2.shape[0]==ang3.shape[0]==a_hmax.shape[0]==a_hmin.shape[0]==a_vert.shape[0])

    mypar['ndir'] = ivdir.shape[0]
    mypar['nst'] = it.shape[0]

    assert (set(it).issubset(set([1,2,3,4,5]))) # is a correct variogram type?

    mypar['ivdir_'] = pd.DataFrame.to_string(pd.DataFrame(ivdir),index= False, header=False) # array to string

    mypar['vst_'] = ''

    for i in range(mypar['ndir']):
        mypar['vst_'] = mypar['vst_'] + str (it[i]) + \
                       '  ' + str (cc[i]) + \
                       '  ' + str (ang1[i]) + \
                       '  ' + str (ang2[i]) + \
                       '  ' + str (ang3[i]) + \
                       '                         - it,cc,ang1,ang2,ang3 \n' +  \
                       '  ' + str (a_hmax[i]) +  \
                       '  ' + str (a_hmin[i]) +  \
                       '  ' + str (a_vert[i]) +  \
                       '                         - a_hmax, a_hmin, a_vert \n'


    par = __vmodel_par.format(**mypar)
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

    # put results in pandas
    nvarg = 1
    ndir = mypar['ndir']
    nlag = mypar['nlag'] + 2

    ignore = np.arange(0,nvarg*ndir*nlag+ndir*nvarg,nlag+1) # list to ignore variogram headers
    # a) read resulting file
    vg = pd.read_csv(mypar['outfl'],
                    header=None,
                    skiprows = ignore,
                    delim_whitespace= True,
                    names = ['Lag',
                            'average separation',
                            'var funct',
                            'number of directions',
                            'covariance ',
                            'correlation'])
    # b) add extra variables from headers
    vg['Variogram'] = np.repeat(range(nvarg), ndir*nlag) # variogram number = row index on parameter['ivpar']
    vg['Direction'] = np.tile(np.repeat(range(ndir), nlag),nvarg)
    vg['tail'] = np.nan
    vg['head'] = np.nan
    vg['type'] = np.nan
    vg['cut'] = np.nan

    # clean a bit zeros and variogram at distance zero
    vg.loc[vg['correlation']==1,'var funct']=None
    vg = vg.set_index(['Variogram', 'Direction', 'Lag'])
    # prepare figure
    fig, ax = plt.subplots(figsize=(8,6))
    for i in vg.index.levels[0]:
        for j in vg.index.levels[1]:
            vg.loc[i,j].plot(kind='line', x= 'average separation', y = 'var funct', ax=ax, label = 'v{} d{}'.format(i,j))

    return vg, fig, ax
