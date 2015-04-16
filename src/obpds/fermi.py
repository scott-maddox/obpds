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

from numpy import exp, log


__all__ = ['fermi_p', 'fermi_n', 'inv_fermi_p', 'inv_fermi_n']


def fermi_p(phi_p, Ev, Nv, Vt):
    return Nv*exp((Ev-phi_p)/Vt)

def fermi_n(phi_n, Ec, Nc, Vt):
    return Nc*exp((phi_n-Ec)/Vt)

def inv_fermi_p(p, Ev, Nv, Vt):
    return Ev - log(p / Nv)*Vt

def inv_fermi_n(n, Ec, Nc, Vt):
    return Ec + log(n / Nc)*Vt