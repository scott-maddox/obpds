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

import numpy


__all__ = ['Solution', 'EquilibriumSolution']


class Solution(object):
    pass

class FlatbandSolution(Solution):
    def __init__(self, T, N, x, Ev, Ec, Ei):
        self.T = T
        self.N = N
        self.x = x
        self.Ev = Ev
        self.Ec = Ec
        self.Ei = Ei

class EquilibriumSolution(Solution):
    def __init__(self, T, N, x, Na, Nd, Ev, Ec, Ei, psi, n, p):
        self.T = T
        self.N = N
        self.x = x
        self.Na = Na
        self.Nd = Nd
        self.Ev = Ev
        self.Ec = Ec
        self.Ei = Ei
        self.psi = psi
        self.n = n
        self.p = p
        self.Ef = numpy.zeros(N)