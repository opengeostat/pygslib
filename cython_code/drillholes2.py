'''
PyGSLIB Drillhole, Module to handle drillhole data and calculations.

Copyright (C) 2015 Adrian Martinez Vargas

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''

# reimplementation of the drillhole module. The old module has some design isues
# that we are planning to fix in this implementation
# This is a pure python module calling optimized.pyx functions for speed

import optimized
