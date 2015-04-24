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

# Make sure we import the local package
import os
import sys
sys.path.insert(0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from obpds import *


from obpds.tests.obpds_test_case import OBPDSTestCase
import unittest

class TestBasic(OBPDSTestCase):
    '''
    Tests basic functionality
    '''

    def test_pn_diode_default(self):
        p = Layer(1*um, Material(GaAs,  1e17/cm3))
        n = Layer(1*um, Material(GaAs, -1e17/cm3))
        d = TwoTerminalDevice(layers=[p, n])
        s = d.get_flatband()
        s = d.get_equilibrium()

    def test_pn_diode_boltzmann(self):
        p = Layer(1*um, Material(GaAs,  1e17/cm3))
        n = Layer(1*um, Material(GaAs, -1e17/cm3))
        d = TwoTerminalDevice(layers=[p, n])
        s = d.get_flatband()
        s = d.get_equilibrium(approx='boltzmann')

    def test_pn_diode_parabolic(self):
        p = Layer(1*um, Material(GaAs,  1e17/cm3))
        n = Layer(1*um, Material(GaAs, -1e17/cm3))
        d = TwoTerminalDevice(layers=[p, n])
        s = d.get_flatband()
        s = d.get_equilibrium(approx='parabolic')

    def test_pn_diode_kane(self):
        p = Layer(1*um, Material(GaAs,  1e17/cm3))
        n = Layer(1*um, Material(GaAs, -1e17/cm3))
        d = TwoTerminalDevice(layers=[p, n])
        s = d.get_flatband()
        s = d.get_equilibrium(approx='kane')

    def test_pn_hj_diode(self):
        p = Layer(1*um, Material(GaAs,  1e17/cm3))
        N = Layer(1*um, Material(AlGaAs(Al=0.3), -1e17/cm3))
        d = TwoTerminalDevice(layers=[p, N])
        s = d.get_flatband()
        s = d.get_equilibrium()

if __name__ == '__main__':
    unittest.main()