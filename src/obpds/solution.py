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
    def __init__(self, T, N, x, Ev, Ec_Gamma, Ec_L, Ec_X, Ec, Ei):
        self.T = T
        self.N = N
        self.x = x
        self.Ev = Ev
        self.Ec_Gamma = Ec_Gamma
        self.Ec_L = Ec_L
        self.Ec_X = Ec_X
        self.Ec = Ec
        self.Ei = Ei

class EquilibriumSolution(Solution):
    def __init__(self, T, N, x, Na, Nd,
                 Ev, Ec_Gamma, Ec_L, Ec_X, Ec, Ei,
                 psi, n_Gamma, n_L, n_X, n, p):
        self.T = T
        self.N = N
        self.x = x
        self.Na = Na
        self.Nd = Nd
        self.Ev = Ev
        self.Ec_Gamma = Ec_Gamma
        self.Ec_L = Ec_L
        self.Ec_X = Ec_X
        self.Ec = Ec
        self.Ei = Ei
        self.psi = psi
        self.n_Gamma = n_Gamma
        self.n_L = n_L
        self.n_X = n_X
        self.n = n
        self.p = p
        self.Ef = numpy.zeros(N)