#
#   Copyright (c) 2015, Scott J Maddox
#
#   This file is part of Open Band Parameters Device Simulator (OBPDS).
#
#   OBPDS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   OBPDS is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with OBPDS.  If not, see <http://www.gnu.org/licenses/>.
#
#############################################################################

import logging; logging.basicConfig()

# Make sure we import the local obpds version
import os
import sys
sys.path.insert(0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from obpds import *

# Layers
n = Layer(1*um, GaAs, -1e16/cm3)

# Device
d = TwoTerminalDevice(layers=[n],
                      contacts=[SchottkyContact(), OhmicContact()],
                      Fn='right')

# Simulate and show the band profile at 0.5 V reverse bias under the zero
# current approximation.
d.show_zero_current(V=-2)

print 'C = {} F/cm**2'.format(d.get_capacitance(V=-0.5))
d.show_cv(-2, 0.2)