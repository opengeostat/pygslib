import pytest
from pygslib.vtktools import inside_polygon, polygon_area


polygon_pt = [ [0., 0., 0.],
               [1., 0., 0.],
               [1., .9, 0.],
               [0., 1., 0.]]

testIn = [[0.5, 0.5, 0.0]]
testOut = [[2.0, 0.5, 0.0]]

def test_inside_polygon_01():
    assert inside_polygon(points = testIn,
                          polygon_pts = polygon_pt,
                          polygon = None,
                          normal = None,
                          bounds = None)[0][0] == 1

    assert inside_polygon(points = testOut,
                          polygon_pts = polygon_pt,
                          polygon = None,
                          normal = None,
                          bounds = None)[0][0] == 0

def test_polygon_area_01():
    assert polygon_area(polygon_pts = polygon_pt,
                        polygon = None) == pytest.approx(0.95, abs = 0.0001)
