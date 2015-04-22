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

# Make sure we import the local obpds version
import os
import sys
sys.path.insert(0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from obpds import *

# Layers
p = Layer(1*um, Material(InAs,  1e15/cm3))
n = Layer(1*um, Material(InAs, -2e18/cm3))

# Device
d = TwoTerminalDevice(layers=[p, n])

# Simulate and show the equilibrium band profile
import matplotlib.pyplot as plt
_, ax1 = plt.subplots()
ax1.set_ymargin(0.05)
ax1.set_ylabel('Energy (eV)')
ax1.set_xlabel('Depth (nm)')

solution = d.get_equilibrium(boltz=True)
x = solution.x*1e7 # nm
ax1.plot(x, solution.Ev, 'r--')
ax1.plot(x, solution.Ec, 'b--', label='Boltzmann')

solution = d.get_equilibrium(boltz=False)
x = solution.x*1e7 # nm
ax1.plot(x, solution.Ev, 'r-')
ax1.plot(x, solution.Ec, 'b-', label='Fermi-Dirac')
ax1.plot(x, solution.Ef, 'k--')

ax1.legend(loc='best')
plt.show()