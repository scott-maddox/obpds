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
from fdint import fd1h, dfd1h, ifd1h


__all__ = ['boltz_p', 'boltz_n', 'iboltz_p', 'iboltz_n',
           'dboltz_p', 'dboltz_n',
           'fermi_p', 'fermi_n', 'ifermi_p', 'ifermi_n',
           'dfermi_p', 'dfermi_n']

# Boltzmann approximation to the Fermi-Dirac integral for a bulk
# semiconductor with parabolic bands.
def boltz_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return exp(phi)*Nv

def boltz_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return exp(phi)*Nc

def iboltz_p(p, Ev, Nv, Vt):
    return Ev - log(p / Nv)*Vt

def iboltz_n(n, Ec, Nc, Vt):
    return Ec + log(n / Nc)*Vt

def dboltz_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return -exp(phi)*Nv/Vt

def dboltz_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return exp(phi)*Nc/Vt

# Fermi-Dirac integral for a bulk semiconductor with parabolic bands
def fermi_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return fd1h(phi)*Nv

def fermi_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return fd1h(phi)*Nc

def ifermi_p(p, Ev, Nv, Vt):
    return Ev - ifd1h(p)*Vt

def ifermi_n(n, Ec, Nc, Vt):
    return Ec + ifd1h(n)*Vt

def dfermi_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return -dfd1h(phi)*Nv/Vt

def dfermi_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return dfd1h(phi)*Nc/Vt
