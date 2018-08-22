# test file to test pytest example: call with ~pytest -v
import pytest

from pygslib.drillhole import ang2cart

#create dum function and assign docstring (pytest does not run directly from ang2cart)
def test_ang2cart_01():
    # use pytest.approx to avoid precision issues
    assert ang2cart(azm = 45, dip = 75) == pytest.approx((0.1830, 0.1830, -0.9659), abs=0.0001)

def test_ang2cart_02():
    assert ang2cart( azm = 0, dip = 0) == pytest.approx((0., 1., 0.), abs = 0.1)

def test_ang2cart_03():
    assert ang2cart( azm = 45, dip = 90) == pytest.approx((0.0, 0.0, -1.0), abs = 0.1)
