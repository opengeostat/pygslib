# -*- coding: utf-8 -*-
#!/usr/bin/env python

# using naming on http://www.gslib.com/gslib_help/programs.html
import subprocess
import copy
import pandas as pd
import pygslib
import numpy as np
import os

__rotcoord_par = \
"""                  Parameters for ROTCOORD
                  ***********************

START OF PARAMETERS:
{datafl}                          -file with data
{icolx} {icoly}                   -columns with X and Y coordinates
{outfl}                           -file for output
{xorigin} {yorigin}               -origin of rotated system in original coordinates (pibot point)
{angle}                           -rotation angle (in degrees clockwise)
{switch}                          -0=convert to rotated coordinate system
                                  -1=convert from rotated system to original system

"""

def rotcoord(parameters, gslib_path = None, silent = False):
    """rotcoord(parameters, gslib_path = None)

    Funtion to rotate coordinates using rotcoord.exe external
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
    pandas.DataFrame with ROTCOORD output

    Example
    --------
    TODO:

    Notes
    ------
    The dictionary with parameters may be as follows::


        parameters = {
            'datafl' : str or, None, or numpy,     # path to file, or none (to use '_xxx_.in') or numpy array (with columns [x,y])
            'icolx'  : int,              # -columns with X and Y coordinates
            'icoly'  : int,
            'outfl'  : str or None,      # path to the output file or None (to use '_xxx_.out')
            'xorigin': float,            # origin of rotated system in original coordinates (pibot point)
            'yorigin': float,
            'angle'  : float,            # rotation angle (in degrees clockwise)
            'switch' : float}            # -0=convert to rotated coordinate system, -1=convert from rotated system to original system


    see http://www.gslib.com/gslib_help/addcoord.html for more information

    """

    if gslib_path is None:
        if os.name == "posix":
            gslib_path = '~/gslib/rotcoord'
        else:
            gslib_path = 'c:\\gslib\\rotcoord.exe'

    mypar = copy.deepcopy(parameters)

    # chek if we use internal files or external and generate files
    if isinstance(parameters['datafl'], np.ndarray):
        # redefine parfile
        mypar['datafl']='_xxx_.in'
        mypar['icolx']=1
        mypar['icoly']=2

        # create dum array
        with open('_xxx_.in',"w") as f:
            f.write('temp file '+'\n')
            f.write('2\n')
            f.write('x\n')
            f.write('y\n')
            np.savetxt(f,parameters['datafl'])
    elif parameters['datafl'] is None:
        mypar['datafl']='_xxx_.in'

    if mypar['outfl'] is None:
        mypar['outfl'] = '_xxx_.out'

    # prepare parameter file and save it
    par = __rotcoord_par.format(**mypar)
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
        raise NameError('gslib rotcoord NameError' + str(stderr.decode('utf-8')))

    if ~silent:
        try:
            print (stdout.decode('utf-8'))
        except:
            print (stdout)

    # return results as panndas array
    return  pygslib.gslib.read_gslib_file(mypar['outfl'])
