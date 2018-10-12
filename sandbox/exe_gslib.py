# -*- coding: utf-8 -*-
#!/usr/bin/env python

# using naming on http://www.gslib.com/gslib_help/programs.html
import subprocess
import copy
import pandas as pd
import pygslib
import numpy as np

__addcoord_par = \
"""                  Parameters for ADDCOORD
                  ***********************

START OF PARAMETERS:
{datafl}                          -file with data
{outfl}                           -file for output
{ireal}                           -realization number
{nx}  {xmn}  {xsiz}               -nx,xmn,xsiz
{ny}  {ymn}  {ysiz}               -ny,ymn,ysiz
{nz}  {zmn}  {zsiz}               -nz,zmn,zsiz
"""

def xstr(s):
    return s or ''

def addcoord(parameters, gslib_path = None):
    """
    Addcoord object.

    Todo: write desc

    Parameters
    ----------
    parameters : dict
        dictionary with parameters
    gslib_path : string (default None)
        absolute or relative path to gslib excecutable programs
    tmp_path : string (default '')
        Path to a temporarty working dir

    Notes
    ------
    The dictionary with parameters may be as follows::


        parameters = {
            'datafl': str,      # path to the input grid file or panda array (full model sorted on IJK required)
            'outfl' : str,      # path to the output grid file or None
            'ireal' : int,      # realization number
            'nx'    : int,      # number of rows, cols and levels
            'ny'    : int,
            'nz'    : int,
            'xmn'   : float,    # coordinates of the centroid of first/corner block
            'ymn'   : float,
            'zmn'   : float,
            'xsiz'  : float,    # grid node separation
            'ysiz'  : float,
            'zsiz'  : float}

    or

            parameters = {
                'datafl': Blockmodel}      # a pygslib Blockmodel


    see http://www.gslib.com/gslib_help/addcoord.html for more information

    """

    if gslib_path is None:
        if os.name == "posix":
            gslib_path = '~/gslib/addcoord'
        else:
            gslib_path = 'c:/gslib/addcoord.exe'

    # chek if we use internal files or external and generate files
    if isinstance(parameters, pygslib.blockmodel.Blockmodel):
        # redefine parfile
        mypar = {}
        mypar['nx'] = parameters.nx
        mypar['ny'] = parameters.ny
        mypar['nz'] = parameters.nz
        mypar['xsiz'] = parameters.dx
        mypar['ysiz'] = parameters.dy
        mypar['zsiz'] = parameters.dz
        mypar['xmn'] = parameters.xorg + parameters.dx/2.
        mypar['ymn'] = parameters.yorg + parameters.dy/2.
        mypar['zmn'] = parameters.zorg + parameters.dz/2.
        mypar['ireal'] = 1

        # check is sorted and full model
        assert ('IJK' in parameters.bmtable.columns)
        assert(all(parameters.bmtable['IJK'].values == np.arange(mypar['nx']*mypar['ny']*mypar['nz'])))

        # update in/out in parfile
        mypar['datafl']='_xxx_.in'
        mypar['outfl'] = None
        # create dum array
        with open('_xxx_.in',"w") as f:
            f.write('temp grid file nx={}, ny={}, nz={}'.format(mypar['nx'],
                                                                mypar['ny'],
                                                                mypar['nz']) +'\n')
            f.write('1\n')
            f.write('IJK\n')
            np.savetxt(f,parameters.bmtable['IJK'].values,fmt='%d')

    else:
        # chek we have all fields
        assert (set(['datafl','outfl','ireal','nx','ny','nz','xmn','ymn','zmn','xsiz','ysiz','zsiz']).issubset(parameters))
        # get a working parameter file
        mypar = copy.deepcopy(parameters)

    if mypar['outfl'] is None:
        mypar['outfl'] = '_xxx_.out'

    # prepare parameter file and save it
    par = __addcoord_par.format(**mypar)
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
        raise NameError('HiThere' + str(stderr.decode('utf-8')))

    print (stdout.decode('utf-8'))

    # return results as panndas array
    return  pygslib.gslib.read_gslib_file(mypar['outfl'])
