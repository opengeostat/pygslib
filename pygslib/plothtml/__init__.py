# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
PyGSLIB.plothtml, Module to plot data in HTML/Jupyter Notebooks.

This is a simplified interface to create javascript plots with
Bokeh and plotly.

The plots implemented are

 - Histograms with GSLIB (suppor declustering weight)
 - CDF with GSLIB (suppor declustering weight)
 - 2D plots with data selection fuctions

Note
----

Copyright 2017, Adrian Martinez Vargas

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE.txt file for details.

"""

# general python modules
import sys
import os.path
import warnings
import pandas as pd
import numpy as np
import bokeh.plotting as bkplt
from bokeh.models import CustomJS, LinearColorMapper, LogColorMapper, ColumnDataSource, ColorBar

# gslib modules
import pygslib.gslib.__plot
import pygslib.gslib.__variograms
import pygslib.gslib.__bigaus
import pygslib.gslib.__bicalib

#pygslib modules
import pygslib

#ipython?
from IPython.display import display, HTML

#-----------------------------------------------------------------------------------------------------------------
#
#    Initialization
#
#-----------------------------------------------------------------------------------------------------------------

# set output for bokeh (Ipython/html file)
if 'ipykernel' in sys.modules:
    bkplt.output_notebook(hide_banner=True)
else:
    bkplt.output_file("tmp_pygslib.html")


#-----------------------------------------------------------------------------------------------------------------
#
#    Plot fuctions (declustering supported)
#
#-----------------------------------------------------------------------------------------------------------------

def show(fig):
    """ Renders a bokeh plot, it just wraps bokeh.io.show()

    """

    bkplt.show(fig)


def histgplt(parameters):
    """Plot gslib.histogram with bokeh, using declustering weight.


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

        parameters = {
            # gslib parameters for histogram calculation
            'hmin' : , # input float (Optional: set as minimum value in dataset). Minimun value in the histogram.
            'hmax' : , # input float (Optional: set as maximum value in dataset). Maximum value in the histogram.
            'ncl'  : , # input int (Optional: set as 10). Number of bins/clases.
            'iwt'  : , # input boolean (Optional: set True). Use weight variable?
            'ilog' : , # input boolean (Optional: set False). If true uses log scale, otherwise uses arithmetic
            'icum' : , # input boolean (Optional: set False). If true uses cumulative histogram, otherwise plots frequency histograms
            'va'   : , # input rank-1 array('d') with bounds (nd). Variable
            'wt'   : , # input rank-1 array('d') with bounds (nd) (Optional, set to array of ones). Declustering weight.
            # visual parameters for figure (if a new figure is created)
            'figure' : , # a bokeh figure object (Optional: new figure created if None). Set none or undefined if creating a new figure.
            'title'  : , # string (Optional, "Histogram"). Figure title
            'xlabel' : , # string (Optional, default "Z"). X axis label
            'ylabel' : , # string (Optional, default "f(%)"). Y axis label
            # visual parameter for the histogram
            'color' : , # string with valid CSS colour (https://www.w3schools.com/colors/colors_names.asp), or an RGB(A) hex value, or tuple of integers (r,g,b), or tuple of (r,g,b,a) (Optional, default "navy")
            'legend': , # string (Optional, default "NA").
            'alpha' : , # float [0-1] (Optional, default 0.5). Transparency of the fill colour
            'lwidth': , # float (Optional, default 1). Line width
            # leyend
            'legendloc': } #  float (Optional, default 'top_right'). Any of top_left, top_center, top_right, center_right, bottom_right, bottom_center, bottom_left, center_left or center


    Returns
    -------
    dict: a dictionary with some data parameters

    {'binval' : , # rank-1 array('d') with bounds (ncl). Class frequency.
     'nincls' : , # rank-1 array('d') with bounds (ncl). Count of values in the class (not using declustering weight).
     'cl'     : , # rank-1 array('d') with bounds (ncl). Class.
     'clwidth': , # rank-1 array('d') with bounds (ncl). Class width
     'xpt025' : , # float. Percentile 0.025
     'xlqt'   : , # float. Percentile 0.25
     'xmed'   : , # float. Mediam (Percentile 0.50)
     'xuqt'   : , # float. Percentile 0.75
     'xpt975' : , # float. Percentile 0.975
     'xmin'   : , # float. Minimum value.
     'xmax'   : , # float. Maximum value.
     'xcvr'   : , # float. Coeficient of variation.
     'xmen'   : , # float. Mean
     'xvar'   : , # float. Variance
     'xfrmx'  : , # float. Maximum Class Frequency
     'dcl'    : } # float. Class width (thmax-thmin)/real(ncl), only usefull if using arithmetic scale.

    fig: a bokeh figure

    Example
    -------

    >>>

    """


    # set weight if not included in the input
    if 'wt' not in parameters:
        parameters['wt']= np.ones(parameters['va'].shape[0], dtype = parameters['va'].dtype)

    if parameters['wt'] is None:
        parameters['wt']= np.ones(parameters['va'].shape[0], dtype = parameters['va'].dtype)

    # set undefine parameters

    if 'ilog' not in parameters:
        parameters['ilog']= 1
    if parameters['ilog'] is None:
        parameters['ilog']= 1
    if parameters['ilog'] < 0:
        parameters['ilog']= 1


    if 'hmin' not in parameters:
        parameters['hmin']= parameters['va'].min()
    if parameters['hmin'] is None:
        parameters['hmin']= parameters['va'].min()
    if parameters['hmin']<=0 and parameters['ilog']==1:
        parameters['hmin']= 0.001

    if 'hmax' not in parameters:
        parameters['hmax']= parameters['va'].max()
    if parameters['hmax'] is None:
        parameters['hmax']= parameters['va'].max()
    if parameters['hmax']<=0 and parameters['ilog']==1:
        raise NameError('All data in parameters["va"]<=0 but you are using logarithmic sale. Make your data positive and try again')

    if 'ncl' not in parameters:
        parameters['ncl']= 10
    if parameters['ncl'] is None:
        parameters['ncl']= 10
    if parameters['ncl'] <= 0:
        parameters['ncl']= 10

    if 'iwt' not in parameters:
        parameters['iwt']= 0
    if parameters['iwt'] is None:
        parameters['iwt']= 0

    if 'icum' not in parameters:
        parameters['icum']= 0
    if parameters['icum'] is None:
        parameters['icum']= 0



    # prepare the output

    parameters_gslib = {'hmin':parameters['hmin'],
                        'hmax':parameters['hmax'],
                        'ncl' :parameters['ncl'],
                        'iwt' :parameters['iwt'],
                        'ilog':parameters['ilog'],
                        'icum':parameters['icum'],
                        'va'  :parameters['va'],
                        'wt'  :parameters['wt']}

    binval,nincls,cl, clwidth,xpt025,xlqt,xmed,xuqt,xpt975, \
    xmin,xmax,xcvr,xmen,xvar,xfrmx,dcl,error = pygslib.gslib.__plot.histplt(**parameters_gslib)

    out1 = {'binval' : binval, # rank-1 array('d') with bounds (ncl). Class frequency.
     'nincls' : nincls , # rank-1 array('d') with bounds (ncl). Count of values in the class (not using declustering weight).
     'cl'     : cl , # rank-1 array('d') with bounds (ncl). Class.
     'clwidth': clwidth, # rank-1 array('d') with bounds (ncl). Class width
     'xpt025' : xpt025, # float. Percentile 0.025
     'xlqt'   : xlqt, # float. Percentile 0.25
     'xmed'   : xmed, # float. Mediam (Percentile 0.50)
     'xuqt'   : xuqt, # float. Percentile 0.75
     'xpt975' : xpt975, # float. Percentile 0.975
     'xmin'   : xmin, # float. Minimum value.
     'xmax'   : xmax, # float. Maximum value.
     'xcvr'   : xcvr, # float. Coeficient of variation.
     'xmen'   : xmen, # float. Mean
     'xvar'   : xvar, # float. Variance
     'xfrmx'  : xfrmx, # float. Maximum Class Frequency
     'dcl'    : dcl} # float. Class width (thmax-thmin)/real(ncl), only usefull if using arithmetic scale.


    # create a figure or update one

    if 'figure' not in parameters:
        parameters['figure']= None

    if 'title' not in parameters:
        if parameters['figure'] is None:
            parameters['title']= 'Histogram'
        else:
            parameters['title']= parameters['figure'].title.text

    if 'xlabel' not in parameters:
        if parameters['figure'] is None:
            parameters['xlabel']= 'Z'
        else:
            parameters['xlabel']= parameters['figure'].xaxis.axis_label

    if 'ylabel' not in parameters:
        if parameters['figure'] is None:
            parameters['ylabel']= 'f(%)'
        else:
            parameters['ylabel']= parameters['figure'].yaxis.axis_label


    if parameters['figure'] is None:
        # create a new plot with a title and axis labels

        if parameters['ilog']==1:
            x_axis_type='log'
        else:
            x_axis_type='linear'

        parameters['figure'] = bkplt.figure(
                x_axis_type=x_axis_type,
                x_range = [min(0.01,parameters['hmin']),  parameters['hmax']], toolbar_location = 'above')

    # set titles
    parameters['figure'].title.text=parameters['title']
    parameters['figure'].xaxis.axis_label=parameters['xlabel']
    parameters['figure'].yaxis.axis_label=parameters['ylabel']


    # add the new plot
    if 'color' not in parameters:
        parameters['color']= 'navy'

    if 'alpha' not in parameters:
        parameters['alpha']= 0.5

    if 'lwidth' not in parameters:
        parameters['lwidth']= 1

    if 'legend' not in parameters:
        parameters['legend']= 'NA'


    parameters['figure'].vbar(x=out1['cl'] -out1['clwidth']/2.,
                              width=-out1['clwidth'], bottom=0,
                              top=out1['binval'],
                              color=parameters['color'],
                              legend = parameters['legend'],
                              alpha=parameters['alpha'],
                              line_width=parameters['lwidth'])

    if 'legendloc' not in parameters:
        parameters['legendloc']= 'top_right'

    parameters['figure'].legend.location = parameters['legendloc']


    #bkplt.show(parameters['figure'])

    return out1, parameters['figure']

def probplt(parameters):
    """Plot gslib.probplt with bokeh, using declustering weight.


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

        parameters = {
            # gslib parameters for histogram calculation
            'iwt'  : , # input boolean (Optional: set True). Use weight variable?
            'va'   : , # input rank-1 array('d') with bounds (nd). Variable
            'wt'   : , # input rank-1 array('d') with bounds (nd) (Optional, set to array of ones). Declustering weight.
            # visual parameters for figure (if a new figure is created)
            'figure' : , # a bokeh figure object (Optional: new figure created if None). Set none or undefined if creating a new figure.
            'title'  : , # string (Optional, "Histogram"). Figure title
            'xlabel' : , # string (Optional, default "Z"). X axis label
            'ylabel' : , # string (Optional, default "f(%)"). Y axis label
            'xlog' : , # boolean (Optional, default True). If true plot X axis in log sale.
            'ylog' : , # boolean (Optional, default True). If true plot Y axis in log sale.
            # visual parameter for the probplt
            'style' : , # string with valid bokeh chart type
            'color' : , # string with valid CSS colour (https://www.w3schools.com/colors/colors_names.asp), or an RGB(A) hex value, or tuple of integers (r,g,b), or tuple of (r,g,b,a) (Optional, default "navy")
            'legend': , # string (Optional, default "NA").
            'alpha' : , # float [0-1] (Optional, default 0.5). Transparency of the fill colour
            'lwidth': , # float (Optional, default 1). Line width
            # leyend
            'legendloc': } #  float (Optional, default 'bottom_right'). Any of top_left, top_center, top_right, center_right, bottom_right, bottom_center, bottom_left, center_left or center


    Returns
    -------
    dict: a dictionary with some data parameters

    {'binval : , # rank-1 array('d') with bounds (nd). Grade class c.
    'cl : , # rank-1 array('d') with bounds (nd). Probability P[z<c] of the grade class c.
    'xpt025' : , # float. Quantile 2.5%
    'xlqt' : , # float . Quantile 25%
    'xmed' : , # float. Median (Quantile 50%)
    'xuqt' : , # float. Quantile 75%
    'xpt975' : , # float. Quantile 97.5%
    'xmin' : , # float. Mininum data value
    'xmax' : , # float. Maximum data value
    'xcvr' : , # float. Coefficient of variation
    'xmen' : , # float. Mean
    'xvar' : , # float. Variance
    'error' : } #int. Error description. If error==0 all is ok.

    fig: a bokeh figure/plot

    Note
    ----
    Common errors are:
    - error = 1  the data has size <= 1
    - error = 10 declustering wight too low

    overlap multiple plots to create special visual effects, e.g. lines with cicles

    Example
    -------

    >>>

    """


    # set weight if not included in the input
    if 'wt' not in parameters:
        parameters['wt']= np.ones(parameters['va'].shape[0], dtype = parameters['va'].dtype)

    if parameters['wt'] is None:
        parameters['wt']= np.ones(parameters['va'].shape[0], dtype = parameters['va'].dtype)

    # set undefine parameters

    if 'iwt' not in parameters:
        parameters['iwt']= 0
    if 'xlog' not in parameters:
        parameters['xlog']= 1
    if 'ylog' not in parameters:
        parameters['ylog']= 1




    # prepare the output

    parameters_gslib = {'iwt':parameters['iwt'],
                        'va':parameters['va'],
                        'wt' :parameters['wt']}

    binval,cl,xpt025,xlqt,xmed,xuqt,xpt975,xmin,xmax, \
    xcvr,xmen,xvar,error = pygslib.gslib.__plot.probplt(**parameters_gslib)

    assert error == 0, 'Error in  pygslib.gslib.__plot.probplt, error = {}'.format('error')

    out1 = {'binval' : binval, # rank-1 array('d') with bounds (nd). Grade class c.
    'cl' : cl, # rank-1 array('d') with bounds (nd). Probability P[z<c] of the grade class c.
    'xpt025': xpt025 , # float. Quantile 2.5%
    'xlqt' : xlqt, # float . Quantile 25%
    'xmed' : xmed, # float. Median (Quantile 50%)
    'xuqt' : xuqt, # float. Quantile 75%
    'xpt975' : xpt975, # float. Quantile 97.5%
    'xmin' : xmin, # float. Mininum data value
    'xmax' : xmax, # float. Maximum data value
    'xcvr' : xcvr, # float. Coefficient of variation
    'xmen' : xmen, # float. Mean
    'xvar' : xvar, # float. Variance
    'error' : error} #int. Error description. If error==0 all is ok.

    # create a figure or update one

    if 'figure' not in parameters:
        parameters['figure']= None

    if 'title' not in parameters:
        if parameters['figure'] is None:
            parameters['title']= 'Probability plot'
        else:
            parameters['title']= parameters['figure'].title.text

    if 'xlabel' not in parameters:
        if parameters['figure'] is None:
            parameters['xlabel']= 'Z'
        else:
            parameters['xlabel']= parameters['figure'].xaxis.axis_label

    if 'ylabel' not in parameters:
        if parameters['figure'] is None:
            parameters['ylabel']= 'Prob [Z<c]'
        else:
            parameters['ylabel']= parameters['figure'].yaxis.axis_label


    if parameters['figure'] is None:
        # create a new plot with a title and axis labels

        if parameters['xlog']:
            x_axis_type='log'
        else:
            x_axis_type='linear'

        if parameters['ylog']:
            y_axis_type='log'
        else:
            y_axis_type='linear'


        parameters['figure'] = bkplt.figure(
                x_axis_type=x_axis_type,
                x_range = [min(0.001, parameters['va'].min()),  parameters['va'].max()], toolbar_location = 'above')

    # set titles
    parameters['figure'].title.text=parameters['title']
    parameters['figure'].xaxis.axis_label=parameters['xlabel']
    parameters['figure'].yaxis.axis_label=parameters['ylabel']


    # add the new plot
    if 'color' not in parameters:
        parameters['color']= 'navy'

    if 'alpha' not in parameters:
        parameters['alpha']= 0.5

    if 'lwidth' not in parameters:
        parameters['lwidth']= 1

    if 'legend' not in parameters:
        parameters['legend']= 'NA'

    # set the plot type

    if 'style' not in parameters:
        parameters['style'] = 'circle'

    if parameters['style'] is None:
        parameters['style'] = 'circle'

    plt = None

    if parameters['style'] == 'circle':
        plt = parameters['figure'].circle
    if parameters['style'] == 'square':
        plt = parameters['figure'].square
    if parameters['style'] == 'cross':
        plt = parameters['figure'].cross
    if parameters['style'] == 'diamond':
        plt = parameters['figure'].diamond
    if parameters['style'] == 'line':
        plt = parameters['figure'].line

    # non of above, then an error
    assert plt is not None, 'parameter style is not valid'

    plt(x=cl, y=binval,
      color=parameters['color'],
      legend = parameters['legend'],
      alpha=parameters['alpha'],
      line_width=parameters['lwidth'])

    if 'legendloc' not in parameters:
        parameters['legendloc']= 'bottom_right'

    parameters['figure'].legend.location = parameters['legendloc']


    #bkplt.show(parameters['figure'])

    return out1, parameters['figure']


def qpplt(parameters):
    """Plot gslib.qpplt with bokeh, using declustering weight.


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

        parameters = {
            # gslib parameters for qq-pp calculation
            'qqorpp': , # integer (Optional, default 0, Q-Q plot). Q-Q plot (qqorpp=0); P-P plot (qqorpp=1)
            'npts'  : , # integer (Optional, default min length of va1 and va2). Number of points to use on the Q-Q or P-P plot (should not exceed the smallest number of data in data1 / data2
            'va1'   : , # rank-1 array('d') with bounds (nd). Variable 1
            'wt1'   : , # rank-1 array('d') with bounds (nd) (Optional, set to array of ones). Declustering weight for variable 1.
            'va2'   : , # rank-1 array('d') with bounds (nd). Variable 2
            'wt2'   : , # rank-1 array('d') with bounds (nd) (Optional, set to array of ones). Declustering weight for variable 2.
            # visual parameters for figure (if a new figure is created)
            'figure' : , # a bokeh figure object (Optional: new figure created if None). Set none or undefined if creating a new figure.
            'title'  : , # string (Optional, "QQ plot" or "PP plot"). Figure title
            'xlabel' : , # string (Optional, default "Z1" or "P1"). X axis label
            'ylabel' : , # string (Optional, default "Z2" or "P2"). Y axis label
            'xlog' : , # boolean (Optional, default True). If true plot X axis in log sale.
            'ylog' : , # boolean (Optional, default True). If true plot Y axis in log sale.
            # visual parameter for the probplt
            'style' : , # string with valid bokeh chart type
            'color' : , # string with valid CSS colour (https://www.w3schools.com/colors/colors_names.asp), or an RGB(A) hex value, or tuple of integers (r,g,b), or tuple of (r,g,b,a) (Optional, default "navy")
            'legend': , # string (Optional, default "NA").
            'alpha' : , # float [0-1] (Optional, default 0.5). Transparency of the fill colour
            'lwidth': , # float (Optional, default 1). Line width
            # leyend
            'legendloc': } #  float (Optional, default 'bottom_right'). Any of top_left, top_center, top_right, center_right, bottom_right, bottom_center, bottom_left, center_left or center


    Returns
    -------
    dict: a dictionary with some data parameters

    {'vr1'  : , # rank-1 array('d') with bounds (npts). Probabilities or Quantiles for variable 1
    'vr2'   : , # rank-1 array('d') with bounds (npts). Probabilities or Quantiles for variable 2
    'error' : } #int. Error description. If error==0 all is ok.

    fig: a bokeh figure/plot

    Note
    ----
    Common errors are:
    - error = 100  ! the number of points may be equal or lowr than the number of points in the smaller dataset
    - error = 10  ! 'Cumulative Probability too LOW'

    Example
    -------

    >>>

    """


    # set weight if not included in the input
    if 'wt1' not in parameters:
        parameters['wt1']= np.ones(parameters['va1'].shape[0], dtype = parameters['va1'].dtype)

    if parameters['wt1'] is None:
        parameters['wt1']= np.ones(parameters['va1'].shape[0], dtype = parameters['va1'].dtype)

    if 'wt2' not in parameters:
        parameters['wt2']= np.ones(parameters['va2'].shape[0], dtype = parameters['va2'].dtype)

    if parameters['wt2'] is None:
        parameters['wt2']= np.ones(parameters['va2'].shape[0], dtype = parameters['va2'].dtype)

    # set npoints
    if 'npoints' not in parameters:
        parameters['npts']= None
    if parameters['npts'] is None:
        parameters['npts']= min(parameters['va1'].shape[0], parameters['va2'].shape[0])

    # set undefine parameters
    if 'qqorpp' not in parameters:
        parameters['qqorpp']= 0   # default QQ plot
    if 'xlog' not in parameters:
        parameters['xlog']= True
    if parameters['xlog'] is None:
        parameters['xlog']= True
    if 'ylog' not in parameters:
        parameters['ylog']= True
    if parameters['ylog'] is None:
        parameters['ylog']= True

    # prepare the output

    parameters_gslib = {'qqorpp':parameters['qqorpp'],
                        'npts':parameters['npts'],
                        'va1':parameters['va1'],
                        'wt1' :parameters['wt1'],
                        'va2':parameters['va2'],
                        'wt2' :parameters['wt2']}

    vr1a,vr2a,error = pygslib.gslib.__plot.qpplt(**parameters_gslib)

    assert error == 0, 'Error in  pygslib.gslib.__plot.qpplt, error = {}'.format('error')

    dmin= min(vr1a.min(),vr2a.min())
    dmax= max(vr1a.max(),vr2a.max())


    out1 = {'vr1a' : vr1a, # rank-1 array('d') with bounds (nd1).
    'vr2a' : vr2a, # rank-1 array('d') with bounds (nd2).
    'error' : error} #int. Error description. If error==0 all is ok.

    # create a figure or update
    if 'figure' not in parameters:
        parameters['figure']= None

    if 'title' not in parameters:
        if parameters['figure'] is None:
            if parameters['qqorpp']==0:
                parameters['title']= 'QQ plot'
            else:
                parameters['title']= 'PP plot'
        else:
            parameters['title']= parameters['figure'].title.text

    if 'xlabel' not in parameters:
        if parameters['figure'] is None:
            if parameters['qqorpp']==0:
                parameters['xlabel']= 'Z1'
            else:
                parameters['xlabel']= 'P[Z1<c]'
        else:
            parameters['xlabel']= parameters['figure'].xaxis.axis_label

    if 'ylabel' not in parameters:
        if parameters['figure'] is None:
            if parameters['qqorpp']==0:
                parameters['ylabel']= 'Z2'
            else:
                parameters['ylabel']= 'P[Z2<c]'
        else:
            parameters['ylabel']= parameters['figure'].yaxis.axis_label


    if parameters['figure'] is None:
        # create a new plot with a title and axis labels

        if parameters['xlog'] is True:
            x_axis_type='log'
        else:
            x_axis_type='linear'
            parameters['xlog']='linear'

        if parameters['ylog'] is True:
            y_axis_type='log'
        else:
            y_axis_type='linear'
            parameters['ylog']='linear'


        parameters['figure'] = bkplt.figure(
                x_axis_type=x_axis_type,
                y_axis_type=y_axis_type,
                x_range = [dmin,  dmax],
                y_range = [dmin,  dmax], toolbar_location = 'above')

    # set titles
    parameters['figure'].title.text=parameters['title']
    parameters['figure'].xaxis.axis_label=parameters['xlabel']
    parameters['figure'].yaxis.axis_label=parameters['ylabel']


    # add the new plot
    if 'color' not in parameters:
        parameters['color']= 'navy'

    if parameters['color'] is None:
        parameters['color']= 'navy'

    if 'alpha' not in parameters:
        parameters['alpha']= 0.5
    if parameters['alpha'] is None:
        parameters['alpha']= 0.5
    if 'lwidth' not in parameters:
        parameters['lwidth']= 1
    if parameters['lwidth'] is None:
        parameters['lwidth']= 1
    if 'legend' not in parameters:
        parameters['legend']= 'NA'

    # set the plot type

    if 'style' not in parameters:
        parameters['style'] = 'circle'
    if parameters['style'] is None:
        parameters['style'] = 'circle'

    plt = None

    if parameters['style'] == 'circle':
        plt = parameters['figure'].circle
    if parameters['style'] == 'square':
        plt = parameters['figure'].square
    if parameters['style'] == 'cross':
        plt = parameters['figure'].cross
    if parameters['style'] == 'diamond':
        plt = parameters['figure'].diamond
    if parameters['style'] == 'line':
        plt = parameters['figure'].line

    #Plot the centerline for reference
    parameters['figure'].line(x=[dmin,dmax], y=[dmin,dmax], color= 'grey', line_width=1, line_dash="4 4")


    # non of above, then an error
    assert plt is not None, 'parameter style is not valid'

    plt(x=vr1a, y=vr2a,
      color=parameters['color'],
      legend = parameters['legend'],
      alpha=parameters['alpha'],
      line_width=parameters['lwidth'])

    if 'legendloc' not in parameters:
        parameters['legendloc']= 'bottom_right'
    if parameters['legendloc'] is None:
        parameters['legendloc']= 'bottom_right'

    parameters['figure'].legend.location = parameters['legendloc']



    return out1, parameters['figure']

#-----------------------------------------------------------------------------------------------------------------
#
#    Plot data with selection supported
#
#-----------------------------------------------------------------------------------------------------------------
def scatter2D(x, y, z, vmin=None, vmax=None, cmax = 'red', cmin = 'blue', cundef='grey',
              cpalete_size=250, ctransform = 'log', radii=1., lcolor= 'black',
              title='Scatterplot',
              TOOLS='resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap,previewsave,box_select,poly_select,lasso_select',
              xlabel='X', ylabel='Y', alpha=0.8, outvar='__inds'):

    """Plot a 2D scaterplot, e.g. data location with bokeh, and allows manual selection of data.


    Parameters
    ----------
        x, y  :  array of floats. Coordinates
        z     :  array of floats. Variable for color
        vmin, vmax : floats (Optional, default None, that is data max and min). Value minimum and maximum for colour scale.
        cmax : string with valid bokeh colour (Optional, default 'red'). Color corresponding to vmax or above
        cmin : string with valid bokeh colour (Optional, default 'blue'). Color corresponding to vmin
        cundef : string with valid bokeh colour (Optional, default 'grey'). Color corresponding to values < vmin and undefined
        cpalete_size : integer (Optional default 250). Size of the colour palete/scale
        ctransform : string ( Optional, default 'log'). If == 'log' will use log scale to asign colors, otherwise will use linear scale
        radii : float (Optional, default 1.). Size of the circles plotted
        lcolor: string (ptional, default 'black'). Color of the circle outline.
        title: string (Optional, default 'Scatterplot'). Title of the plot
        TOOLS: string (Optional, default 'resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap,previewsave,box_select,poly_select,lasso_select'). Tools shown in the plot.
        xlabel: string (Optional, default 'X'). X axis label
        ylabel: string (Optional, default 'Y'). Y axis label
        alpha: float (Optional, default 0.8). Transparency of circles background
        outvar: string (Optional, default '__inds'). Name of the callback variable. Each time you select points in the scatterplot this variable will be updated.


    Note
    -------
    This function will create a new variable defined by outvar, e.g. '__inds', the firts time you select points. This variable will be overwritten if it already exist.
    This variable may be overwritten by selecting points in a different plot with same value in the parameter outvar.

    The selection variable is updated in your python kernel from a javascript (web brouser) using a calling back function.


    Example
    -------
    >>> pygslib.plothtml.scatter2D(
              # data
              x= mydata['Xlocation'],
              y= mydata['Ylocation'],
              z= mydata['Primary'],
              # vmin and max for formating
              vmin=0.1,
              vmax=30,
              # parameters for color
              cmax = 'red',
              cmin = 'blue',
              cundef='grey',
              cpalete_size=300,
              ctransform = 'linear',
              #parameter for circle size
              radii=0.5,
              # line color
              lcolor= 'black',
              # parameters for plot
              title='Scatterplot',
              TOOLS='resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap,previewsave,box_select,poly_select,lasso_select',
              xlabel='X',
              ylabel='Y',
              alpha=0.9,
              # call back variable
              outvar='__inds')
    >>> print " Select some datapoints, e.g. using the tool lasso "
    >>> print mydata.iloc[__inds,:]

    """

    #prepare data
    df = pd.DataFrame({'x': x, 'y':y, 'z':z})
    df['radii']=radii
    df['alpha']=alpha
    df['lcolor'] = lcolor

    if vmin is None:
        vmin = z.min()
    if vmax is None:
        vmax = z.max()

    assert vmin<vmax

    source = ColumnDataSource(df)


    # create a palete.
    colourscale=pygslib.charttable.Leyend_num(vmin=vmin, vmax=vmax, cmax = cmax, cmin = cmin,
                                              undef = cundef, nval= cpalete_size, convert='HEX')
    palete = colourscale.c.tolist()

    # create a mapper
    if ctransform == 'log':
       assert vmin > 0 and vmax > 0, "vmin/vmax <=0 and using log transform"
       mapper = LogColorMapper(palette=palete, low=vmin, high=vmax)
    else:
       mapper = LinearColorMapper(palette=palete, low=vmin, high=vmax)

    mapper.low_color = cundef
    mapper.high_color = cmax

    # create figure
    p = bkplt.figure(title = title, tools=TOOLS, toolbar_location = 'above')
    p.xaxis.axis_label = xlabel
    p.yaxis.axis_label = ylabel

    # colorbar
    color_bar = ColorBar(color_mapper=mapper, location=(0, 0))

    # plot
    p.scatter(source=source,
              x="x", y="y",radius='radii',
              fill_color={'field': 'z', 'transform': mapper},
              fill_alpha='alpha',
              line_color='lcolor')

    p.add_layout(color_bar, 'right')


    source.callback = CustomJS(args=dict(p=p), code="""
            var inds = cb_obj.get('selected')['1d'].indices;
            var d1 = cb_obj.get('data');
            console.log(d1)
            var kernel = IPython.notebook.kernel;
            IPython.notebook.kernel.execute("{} = " + inds);
            """.format(outvar)
    )

    bkplt.show(p)
