"""
Charttable: HTML Tables with bars charts for exploratory data analysis

Copyright 2016, Adrian Martinez Vargas.
Licensed under GPL.
"""


"""
Some notes
The HTML code of this module was inspired on 
http://alistapart.com/article/accessibledatavisualization#horizontal_bar_charts 

The css for extra formating is not included, Using the global parameters from IPYTHON

"""

import warnings
from IPython.display import display, HTML 
from colour import Color
import numpy as np

class Leyend_num(object):
    """
    Leyend_num(vmin=0, 
               vmax=10, 
               cmax = Color("red"), 
               cmin = Color("blue"), 
               undef = Color("grey"), 
               nval= 5)
    
    Define a legend class for continuous numerical intervals
    
    Parameters
    ----------
    vmin, vmax : numerical value 
        Minimum and maximum values for the colour scale. 
    cmin, cmax : instances of class Color  
        Minimum and maximum colours for the colour scale. 
    undef : instance of class Color 
        Color for undefined values (nan or our of range). 
    nval : numerical value 
        Number of intervals in the colour scale. 
           
    Attributes
    ----------
    c : Numpy array of Color instances 
        colours of the colour scale. c[0] is for undefined values
    v : Numpy array of numeric values 
        corresponding value interval [v(i),v(i+1)[ for the colour c(i).
        [v(n), inf[ will take the colour c(n). 
        
    Note
    ----
    
    Use function apply_colour to get the color for each value in an 
    arbitrary array
    
    """    
    # general parameteres
    nval= None
    vmin= None
    vmax= None
    cmax= None 
    cmin= None
    undef= None
    
    # value array
    c= None
    v= None
    
    def __init__(self, vmin=0, vmax=10, cmax = "red", cmin = "blue", undef = "grey", nval= 5, convert=None, log=True): 

        #check parameters
        #assert isinstance(cmax, Color), " Error : cmax is not an instance of Color "
        #assert isinstance(cmin, Color), " Error : cmin is not an instance of Color "
        #assert isinstance(undef, Color), " Error : undef is not an instance of Color "

        # general parameters
        self.nval= nval
        self.vmin=vmin
        self.vmax=vmax
        self.cmax = Color(cmax) 
        self.cmin = Color(cmin)
        self.undef = Color(undef)

        # define arrays
        self.c= np.empty(self.nval+1, dtype = object)
        self.v= np.empty(self.nval+1, dtype = object) 
        
        # set undefined colour        
        self.c[0]= self.undef
        self.v[0]= np.nan   
        
        # asign values to color and value interval
        self.c[1:]= np.array(list(self.cmin.range_to(self.cmax, self.nval)))
        
        if log is True:
            self.v[1:]= np.logspace(np.log(self.vmin), np.log(self.vmax), self.nval, base = np.e)
        else:    
            self.v[1:]= np.linspace(self.vmin, self.vmax, self.nval)
            
        if convert=='HEX':
            for i in range(self.c.shape[0]):
                self.c[i]=self.c[i].hex_l
    
    
    def apply_colour(self, data): 
        """
        Returns colour corresponding to each value in an array
        
                
        Examples
        --------
        
        >>> cc = my_leyend_num.apply_colour([77, 87.9, 77])
         
        """
             
        cc = np.empty(data.shape[0], dtype = object)
        cc[:] = np.nan
        
        # assign the nan
        mask = (data < self.v.min()) | (~np.isfinite(data))
        cc[mask] = self.c[0]
        
        # assign other colours
        for i in range(1, self.nval):
            mask = (data >= self.v[i]) & (data < self.v[i+1]) 
            cc[mask] = self.c[i]
        
        # las value
        mask = (data >= self.v[self.nval])
        cc[mask] = self.c[self.nval]
        
        return cc


class Leyend_cat(object):
    """
    Leyend_cat(values 
               cmax = Color("red"), 
               cmin = Color("blue"), 
               undef = Color("grey"), 
               nval= 5)
    
    Define a legend class for categories
    
    Parameters
    ----------
    values : array  
        Categories with any arbitrary value. 
    cmin, cmax : instances of class Color  
        Minimum and maximum colours for the colour scale. 
    undef : instance of class Color 
        Color for undefined values (nan or our of range). 
    nval : numerical value 
        Number of intervals in the colour scale. 
           
    Attributes
    ----------
    c : Numpy array of Color instances 
        colours of the colour scale. c[0] is for undefined values
    v : Numpy array of numeric values 
        corresponding value interval [v(i),v(i+1)[ for the colour c(i).
        [v(n), inf[ will take the colour c(n). 
        
    Note
    ----
    
    Use function apply_colour to get the color for each value in an 
    arbitrary array
    
    """        
    # general parameteres
    nval= None
    cmax= None 
    cmin= None
    undef= None
    
    # value array
    c= None
    v= None
    
    def __init__(self, values, cmax = "red", cmin = "blue", undef = "grey", convert=None): 

        #check parameters
        #assert isinstance(cmax, Color), " Error : cmax is not an instance of Color "
        #assert isinstance(cmin, Color), " Error : cmin is not an instance of Color "
        #assert isinstance(undef, Color), " Error : undef is not an instance of Color "
        
        # general parameters
        self.nval = len(values)
        self.cmax = Color(cmax)
        self.cmin = Color(cmin)
        self.undef = Color(undef)

        # define arrays
        self.v = np.empty(self.nval+1, dtype = object)
        self.c = np.empty(self.nval+1, dtype = object)
        
        # set undefined colour  
        self.v[0] = np.nan
        self.c[0] = Color(undef)
        
        # asign values to color and value interval
        self.v[1:] = np.sort(values)
        self.c[1:] = np.array(list(self.cmin.range_to(self.cmax, self.nval)))
        
        if convert=='HEX':
            for i in range(self.c.shape[0]):
                self.c[i]=self.c[i].hex_l


    def apply_colour(self, data): 
        """
        Returns colour corresponding to each value in an array
        
                
        Examples
        --------
        
        >>> cc = my_leyend_num.apply_colour(['granite', 'basalt', 'NA'])
         
        """        
        dat = np.array(data)
        
        cc = np.empty(data.shape[0], dtype = object)
        cc[:] = np.nan
        
        # assign the nan
        cc[:] = self.c[0]
        
        # assign other colours
        for i in range(1, self.nval+1): 
            cc[dat==self.v[i]] = self.c[i]
            
        
        return cc
 
 

class HTMLTable(object):
    """
    HTMLTable( data )
    
    Define and create tables with special formats: bars and coloured
    columns. This is useful for exploratory data analysis
   
           
    Attributes
    ----------
    table : dict 
        Table with special columns
    html : str
        HTML code of the table with special columns
    
    nrows : int
        Number of rows
    
    columns: list of str
        Columns labels. The table will be plot in the same order 
    
    Note
    ----
    
    This is intended to plot nice tables with special format in Ipython
    or Jupyter notebooks. The formats are
    
    a) simple text
    b) colored cells
    c) colored bars   
    
    The parameters for colour, bar size and text can be different, for
    example, a color bar can contain color as a rock type, bar size as 
    gold grade and text as Silver grade. In other words 3 different source
    of information can be represented in one column
    
    """      
    
    
    # css for bars configuration
    __css = """

      <style type="text/css">
        table {
          border-collapse: collapse;
        }

        table, th, td {
          border: 1px solid black;
        }
      
        /* CHART LISTS */
        .chartlist { 
          float: left; 
          border-top: 0px solid #EEE; 
          width: 100%;
        }
        .chartlist .chart { 
          position: relative;
          margin: 0 0 0 0;
          border-bottom: 0px solid #EEE; 
        }
        .chartlist .chart p { 
          display: block; 
          padding: 0.0em 0em 0.0em 0.0em;
          position: relative; 
          z-index: 2; 
        }
        .chartlist .index { 
          display: block; 
          position: absolute; 
          top: 0; 
          left: 0; 
          height: 100%; 
          text-indent: -9999px; 
          overflow: hidden; 
          line-height: 1em;
        }
      </style>
    """    
    
    
    def __init__ (self):
        #properties   
        self.table = {}     
        self.nrows = 0 

        self.columns = []

        self.html = ""
    
    def clear_table(self):
        """
        Removes all data from self.Table and reassign index from self.data 
        
        """
        #clear table
        self.table={}
                                
        self.nrows = 0
     
         
                     
    
    def addcolumn(self, name, text=None, bar=None, colour=None, overwrite=False): 
        """
        addcolumn(name, text, bar=None, colour=None,overwrite=False)
        
        Add a special column to self.table
        
        Parameters
        ----------
        name : str 
            column name 
        text : 1D array like
            array of textual values for the column
        bar : 1D array like 
            array with bar sizes (0-100) 
        colour : 1D array like 
            array of colours compatibles to the class Colour
        overwrite: boolean
            set to True to overwrite column with existing name
        
        """
        
        
        
        #Check that we have data
        assert not all([(text is None), (bar is None), (colour is None)]), 'you may provide data' 
                
        # check that the new column name is not in the table 
        if overwrite==False:
            assert name not in self.table.keys(), 'The name "{}" already exist in table, use overwrite=True to replace'.format(name) 

        cformat = [False,False,False]
        
        #check some format and data
        if text is not None: 
            cformat[0] = True 
            if self.nrows>0: # adding new column
                assert len(text) == self.nrows, 'Bad shape in text len(text) != {} '.format(self.nrows)
            else:
                self.nrows= len(text)

            
        if bar is not None: 
            cformat[1] = True 
            if self.nrows>0:
                assert len(bar)  == self.nrows, 'Bad shape in text len(text) != {} '.format(self.nrows)
            else:
                self.nrows= len(bar)            
                if text is not None:
                    assert len(text) == len(bar), 'Bad shape len(text) != len(bar) '
            
            
        if colour is not None: 
            cformat[2] = True
            ncol=len(text)
            if self.nrows>0:
                assert len(colour)==self.nrows, 'Bad shape in colour; arrays with len != {} '.format(self.nrows)
            else:
                self.nrows= len(colour)
                if text is not None:
                    assert len(text) == len(colour), 'Bad shape len(text) != len(colour) '
                if bar is not None:
                    assert len(bar) == len(colour), 'Bad shape len(bar) != len(colour) '          
        
        #create the template dictionary
        column = {'name': name,  # name to be displayed in table HTML
                  'format': cformat , # format [Text(t/f), bar(t/f), ccolour(t/f)]
                  'text':text, # text to be diplayed in each cell
                  'bar': bar,                   # bar size [0-100%] to be plotted
                  'colour': colour}               # colour asigned to the bar
        
        
        self.table[name]= column
        
        self.columns.append(name)
        
    
    def make_html(self):
        """
        make_html()
        
        Makes the HTML code of the table on self.html
        
        """
        
        self.html = ""
        
        # add header
        self.html = self.__css + "<table class='chartlist'>\n"
        
        #add column header
        self.html += "    <tr>\n"
        for c in self.table.keys():
            self.html += "        <th>{0}</td>\n".format(c)
        
        #close header row
        self.html += "    <tr>\n"
        
        
        # add data row by row
        for i in range(self.nrows):
            
            
            # add row
            self.html += "    <tr>\n" 
            
            #add columns
            for c in self.columns:

                # get data
                text = ""
                if self.table[c]['format'][0]:
                    text = self.table[c]['text'][i]
                
                backgr = ""
                if self.table[c]['format'][2]: 
                    colour= Color(self.table[c]['colour'][i])                    
                    backgr = " background: {0};".format(colour.hex_l)
                
                bar = "width: 100%;"
                if self.table[c]['format'][1]:
                    br = self.table[c]['bar'][i] 
                    if br < 0: br=0
                    if br >100: br=100
                    
                    bar = "width: {0}%;".format (br)
            
    
                self.html += "      <td class='chart'>\n"            
                self.html += "            <p>{0}</p>\n".format(text)
                self.html += "            <span class='index' style=' {0} {1}'></span>\n".format(bar, backgr)
                self.html += "      </td>\n"
                    
                    
            
            #close row
            self.html += "    </tr>\n"
            
        #close table
        self.html += "</table>\n"
        
    def display(self):
        """
        display()
        
        Display the Table in Ipython/Jupyter notebooks
        
        """        
        
        self.make_html()
        
        display(HTML(self.html))

