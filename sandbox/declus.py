# -*- coding: utf-8 -*-
#!/usr/bin/env python

# using naming on http://www.gslib.com/gslib_help/programs.html
import subprocess
import copy
import pandas as pd
import pygslib
import numpy as np
import os

__declus_par = \
"""                  Parameters for DECLUS
                  *********************

START OF PARAMETERS:
{datafl}                          -file with data
{icolx} {icoly} {icolz} {icolvr}                   -columns for X, Y, Z and variable
{tmin} {tmax}               -  trimming limits
{sumfl}                     -file for summary output
{outfl}                     -file for output with data & weights
{anisy} {anisz}             -Y and Z cell anisotropy (Ysize=size*Yanis)
{minmax}                    -0=look for minimum declustered mean (1=max)
{ncell} {cmin} {cmax}       -number of cell sizes, min size, max size
{noff}                      -number of origin offsets
"""

def declus(parameters, gslib_path = None, silent = False):
    """declus(parameters, gslib_path = None)

    Funtion to do declustering using declus.exe external
    gslib program.

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
    (pandas.DataFrame, pandas.DataFrame) with declus output and declus summary

    Example
    --------
    TODO:

    Notes
    ------
    The dictionary with parameters may be as follows::


        parameters = {
            'datafl' : str or, None, or numpy,     # path to file, or none (to use '_xxx_.in') or numpy array (with columns [x,y])
            'icolx'  : int,                        # -columns for X, Y, Z and variable
            'icoly'  : int,
            'icolz'  : int,
            'icolvr' : int,
            'tmin'   : float,            # trimming limits min and max (raws out of this range will be ignored)
            'tmax'   : float,
            'sumfl'  : str or None,      # path to the output summary file or None (to use '_xxs_.out')
            'outfl'  : str or None,      # path to the output file or None (to use '_xxx_.out')
            'anisy': float,              # Y and Z cell anisotropies (Ysize=size*Yanis)
            'anisz': float,
            'minmax' : int,              # 0=look for minimum declustered mean (1=max)
            'ncell' : int,               # number of cell sizes, min size, max size
            'cmin' : float,
            'cmax' : float,
            'noff' : int}                # number of origin offsets


    see http://www.gslib.com/gslib_help/declus.html for more information

    """

    if gslib_path is None:
        if os.name == "posix":
            gslib_path = '~/gslib/declus'
        else:
            gslib_path = 'c:\\gslib\\declus.exe'

    mypar = copy.deepcopy(parameters)

    # chek if we use internal files or external and generate files
    if isinstance(parameters['datafl'], np.ndarray):
        # redefine parfile
        mypar['datafl']='_xxx_.in'
        mypar['icolx']=1
        mypar['icoly']=2
        mypar['icolz']=3
        mypar['icolvr']=4

        # create dum array
        with open('_xxx_.in',"w") as f:
            f.write('temp file '+'\n')
            f.write('4\n')
            f.write('x\n')
            f.write('y\n')
            f.write('z\n')
            f.write('variable\n')
            np.savetxt(f,parameters['datafl'])
    elif parameters['datafl'] is None:
        mypar['datafl']='_xxx_.in'

    if mypar['outfl'] is None:
        mypar['outfl'] = '_xxx_.out'

    if mypar['sumfl'] is None:
        mypar['sumfl'] = '_xxs_.out'

    # prepare parameter file and save it
    par = __declus_par.format(**mypar)
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
    return  pygslib.gslib.read_gslib_file(mypar['outfl']), pygslib.gslib.read_gslib_file(mypar['sumfl'])
