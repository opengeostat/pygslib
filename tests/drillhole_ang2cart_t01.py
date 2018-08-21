# test file: call with ~pytest --doctest-modules
#Load function and its docstring
from pygslib.drillhole import ang2cart
from pygslib.drillhole import Drillhole
#import doctest
#doctest.run_docstring_examples(ang2cart,locals())

#create dum function and assign docstring (pytest does not run directly from ang2cart)
def dum():
    pass

dum.__doc__ = ang2cart.__doc__
