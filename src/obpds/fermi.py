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
from numpy import exp, log, sqrt, pi
from fdint import (parabolic, dparabolic, iparabolic,
                   nonparabolic, dnonparabolic,
                   vparabolic, vdparabolic, viparabolic,
                   vnonparabolic, vdnonparabolic)


__all__ = ['boltz_p', 'boltz_n',
           'iboltz_p', 'iboltz_n',
           'boltz_dp', 'boltz_dn',
           'parabolic_p', 'parabolic_n',
           'iparabolic_p', 'iparabolic_n',
           'parabolic_dp', 'parabolic_dn',
           'nonparabolic_n', 'nonparabolic_dn']

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

def boltz_dp(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return -exp(phi)*Nv/Vt

def boltz_dn(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return exp(phi)*Nc/Vt

# Fermi-Dirac integral for a bulk semiconductor with parabolic bands
assert 2/sqrt(pi) == 1.1283791670955126
def parabolic_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return parabolic(phi)*(Nv*1.1283791670955126)

def parabolic_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return parabolic(phi)*(Nc*1.1283791670955126)

def iparabolic_p(p, Ev, Nv, Vt):
    return Ev - iparabolic(p / (Nv*1.1283791670955126))*Vt

def iparabolic_n(n, Ec, Nc, Vt):
    return Ec + iparabolic(n / (Nc*1.1283791670955126))*Vt

def parabolic_dp(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return dparabolic(phi)*(Nv*-1.1283791670955126)/Vt

def parabolic_dn(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return dparabolic(phi)*(Nc*1.1283791670955126)/Vt
    

# Fermi-Dirac integral for a bulk semiconductor with nonparabolic bands
def nonparabolic_n(phi_n, Ec, Nc, alpha, Vt):
    '''
    Nonparabolic Fermi-Dirac integral.
    '''
    phi = (phi_n-Ec)/Vt
    return nonparabolic(phi, alpha)*(Nc*1.1283791670955126)

def nonparabolic_dn(phi_n, Ec, Nc, alpha, Vt):
    '''
    Derivative of the nonparabolic Fermi-Dirac integral.
    '''
    phi = (phi_n-Ec)/Vt
    return dnonparabolic(phi, alpha)*(Nc*1.1283791670955126)/Vt

if __name__ == "__main__":
    from scipy.integrate import quad
    def _num_fermi(phi, alpha):
        result = quad(lambda x: sqrt(x)/(1.+exp(x-phi)),
                      0., numpy.inf)[0]
        return result
    num_fermi = numpy.vectorize(_num_fermi)
    def _num_npfermi(phi, alpha):
        result = quad(lambda x: sqrt((x)*(1.+alpha*x))*
                                (1.+2.*alpha*x)/(1.+exp(x-phi)),
                      0., numpy.inf)[0]
        return result
    num_npfermi = numpy.vectorize(_num_npfermi)

    import matplotlib.pyplot as plt
    _, ax = plt.subplots()
    # boltzmann
    phi = numpy.linspace(-30, 10, 10000)
    ax.semilogy(phi, exp(phi)/1.1283791670955126, 'b')

    # parabolic
    phi = numpy.linspace(-30, 100, 1000)
    ax.semilogy(phi, parabolic(phi), 'r')
    def _num_fd1h(phi):
        result = quad(lambda x: sqrt(x)/(1.+exp(x-phi)), 0., 100.)[0]
        return result
    num_fd1h = numpy.vectorize(_num_fd1h)
    ax.semilogy(phi, num_fd1h(phi), 'r--', lw=2)

    # non-parabolic
    phi = numpy.linspace(-30, 100, 1000)
    alpha = numpy.empty_like(phi)
    alpha.fill(0.07)
    ax.semilogy(phi, nonparabolic(phi, alpha), 'g')
    ax.semilogy(phi, dnonparabolic(phi, alpha), 'g:')
    ax.semilogy(phi, num_npfermi(phi, alpha), 'g--', lw=2)

    plt.show()