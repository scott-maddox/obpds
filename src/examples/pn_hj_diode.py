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
from obpds import GaAs, AlGaAs, cm3, um, Material, Layer, OhmicContact, LayerStructure

# Layers
p = Layer(1*um, Material(GaAs,  1e17/cm3))
N = Layer(1*um, Material(AlGaAs(Al=0.3), -1e17/cm3))

# Contacts
top = OhmicContact()
bottom = OhmicContact()

# Layer Structure
ls = LayerStructure([top, p, N, bottom])

# ls.show_composition() # show the composition vs. depth
# ls.show_doping() # show the doping vs. depth
# ls.show_flatband() # show the flatband profile vs. depth

# Simulate and show the equilibrium band profile using the default method.
ls.show_equilibrium(N=1000)