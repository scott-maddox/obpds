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
from .fd import para, para_d, ipara, kane, kane_d


__all__ = ['boltz_p', 'boltz_n', 'boltz_phi_p', 'boltz_phi_n',
           'boltz_dp', 'boltz_dn',
           'para_p', 'para_n', 'ipara_p', 'ipara_n',
           'para_dp', 'para_dn',
           'kane_n', 'kane_dn']

# Boltzmann approximation to the Fermi-Dirac integral for a bulk
# semiconductor with parabolic bands.
def boltz_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return exp(phi)*Nv

def boltz_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return exp(phi)*Nc

def boltz_phi_p(p, Ev, Nv, Vt):
    return Ev - log(p / Nv)*Vt

def boltz_phi_n(n, Ec, Nc, Vt):
    return Ec + log(n / Nc)*Vt

def boltz_dp(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return -exp(phi)*Nv/Vt

def boltz_dn(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return exp(phi)*Nc/Vt

# Fermi-Dirac integral for a bulk semiconductor with parabolic bands
assert 2/sqrt(pi) == 1.1283791670955126
def para_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return para(phi)*(Nv*1.1283791670955126)

def para_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return para(phi)*(Nc*1.1283791670955126)

def ipara_p(p, Ev, Nv, Vt):
    return Ev - ipara(p / (Nv*1.1283791670955126))*Vt

def ipara_n(n, Ec, Nc, Vt):
    return Ec + ipara(n / (Nc*1.1283791670955126))*Vt

def para_dp(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return para_d(phi)*(Nv*-1.1283791670955126)/Vt

def para_dn(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return para_d(phi)*(Nc*1.1283791670955126)/Vt
    
def kane_n(phi_n, Ec, Nc, alpha, Vt):
    '''
    Non-parabolic Fermi-Dirac integral.
    '''
    phi = (phi_n-Ec)/Vt
    return kane(phi, alpha)*(Nc*1.1283791670955126)

def kane_dn(phi_n, Ec, Nc, alpha, Vt):
    '''
    Derivative of the non-parabolic Fermi-Dirac integral.
    '''
    phi = (phi_n-Ec)/Vt
    return kane_d(phi, alpha)*(Nc*1.1283791670955126)/Vt

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
    ax.semilogy(phi, para(phi), 'r')
    def _num_fd1h(phi):
        result = quad(lambda x: sqrt(x)/(1.+exp(x-phi)), 0., 100.)[0]
        return result
    num_fd1h = numpy.vectorize(_num_fd1h)
    ax.semilogy(phi, num_fd1h(phi), 'r--', lw=2)

    # non-parabolic
    phi = numpy.linspace(-30, 100, 1000)
    alpha = numpy.empty_like(phi)
    alpha.fill(0.07)
    ax.semilogy(phi, kane(phi, alpha), 'g')
    ax.semilogy(phi, kane_d(phi, alpha), 'g:')
    ax.semilogy(phi, num_npfermi(phi, alpha), 'g--', lw=2)

    plt.show()