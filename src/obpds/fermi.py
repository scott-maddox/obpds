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
from fdint import fdk, dfdk
from ifdint import ifd1h


__all__ = ['boltz_p', 'boltz_n', 'iboltz_p', 'iboltz_n',
           'dboltz_p', 'dboltz_n',
           'fermi_p', 'fermi_n', 'ifermi_p', 'ifermi_n',
           'dfermi_p', 'dfermi_n',
           'npfermi', 'npfermi_n', 'dnpfermi_n']

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
assert 2/sqrt(pi) == 1.1283791670955126
def fermi_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return fdk(0.5,phi)*(Nv*1.1283791670955126)

def fermi_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return fdk(0.5,phi)*(Nc*1.1283791670955126)

def ifermi_p(p, Ev, Nv, Vt):
    return Ev - ifd1h(p / (Nv*1.1283791670955126))*Vt

def ifermi_n(n, Ec, Nc, Vt):
    return Ec + ifd1h(n / (Nc*1.1283791670955126))*Vt

def dfermi_p(phi_p, Ev, Nv, Vt):
    phi = (Ev-phi_p)/Vt
    return -0.5*fdk(-0.5,phi)*(Nv*1.1283791670955126)/Vt

def dfermi_n(phi_n, Ec, Nc, Vt):
    phi = (phi_n-Ec)/Vt
    return 0.5*fdk(-0.5,phi)*(Nc*1.1283791670955126)/Vt

def _npfermi(phi, alpha):
    '''
    Approximation of the Fermi-Dirac integral for a bulk semiconductor with
    a non-parabolic band. This approximation degrades significantly for
    alpha > 0.07, particularly at phi ~= 15.
    '''
    if phi < 20:
        # taylor series approximation around alpha = 0
        return (fdk(0.5,phi)+
                 alpha*(2.5*fdk(1.5,phi)+
                 alpha*(0.875*fdk(2.5,phi)+
                 alpha*(-0.1875*fdk(3.5,phi)+
                 alpha*(0.0859375*fdk(4.5,phi)+
                 alpha*(-0.05078125*fdk(5.5,phi)+
                 alpha*(0.0341796875*fdk(6.5,phi)+
                 alpha*(-0.02490234375*fdk(7.5,phi)+
                 alpha*(0.019134521484375*fdk(8.5,phi)+
                 alpha*(-0.0152740478515625*fdk(9.5,phi)+
                 alpha*(0.012546539306640625*fdk(10.5,phi)
                        )))))))))))
    else:
        # sommerfeld approximation for phi -> inf
        return 0.6666666666666666*(phi*(1.+alpha*phi))**1.5
npfermi = numpy.vectorize(_npfermi)

def _dnpfermi(phi, alpha):
    '''
    Approximation of the derivative of the Fermi-Dirac integral for a bulk
    semiconductor with a non-parabolic band. This approximation degrades
    significantly for alpha > 0.07, particularly at phi ~= 15.
    '''
    if phi < 20:
        # taylor series approximation around alpha = 0
        return (dfdk(0.5,phi)+
                 alpha*(2.5*dfdk(1.5,phi)+
                 alpha*(0.875*dfdk(2.5,phi)+
                 alpha*(-0.1875*dfdk(3.5,phi)+
                 alpha*(0.0859375*dfdk(4.5,phi)+
                 alpha*(-0.05078125*dfdk(5.5,phi)+
                 alpha*(0.0341796875*dfdk(6.5,phi)+
                 alpha*(-0.02490234375*dfdk(7.5,phi)+
                 alpha*(0.019134521484375*dfdk(8.5,phi)+
                 alpha*(-0.0152740478515625*dfdk(9.5,phi)+
                 alpha*(0.012546539306640625*dfdk(10.5,phi)
                        )))))))))))
    else:
        # sommerfeld approximation for phi -> inf
        return (2. * (phi*(alpha*phi + 1))**1.5 * 
                (1.0*alpha*phi + 0.5) / (phi*(alpha*phi + 1)))
dnpfermi = numpy.vectorize(_dnpfermi)
    
def npfermi_n(phi_n, Ec, Nc, alpha, Vt):
    '''
    Non-parabolic Fermi-Dirac integral.
    '''
    phi = (phi_n-Ec)/Vt
    return npfermi(phi, alpha)*(Nc*1.1283791670955126)

def dnpfermi_n(phi_n, Ec, Nc, alpha, Vt):
    '''
    Derivative of the non-parabolic Fermi-Dirac integral.
    '''
    phi = (phi_n-Ec)/Vt
    return dnpfermi(phi, alpha)*(Nc*1.1283791670955126)/Vt

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

    import numpy
    import matplotlib.pyplot as plt
    _, ax = plt.subplots()
    # boltzmann
    phi = numpy.linspace(-30, 10, 10000)
    ax.semilogy(phi, exp(phi)/1.1283791670955126, 'b')

    # parabolic
    phi = numpy.linspace(-30, 100, 1000)
    ax.semilogy(phi, fdk(0.5,phi), 'r')
    from scipy.integrate import quad
    def _num_fd1h(phi):
        result = quad(lambda x: sqrt(x)/(1.+exp(x-phi)), 0., 100.)[0]
        return result
    num_fd1h = numpy.vectorize(_num_fd1h)
    ax.semilogy(phi, num_fd1h(phi), 'r--', lw=2)

    # non-parabolic
    phi = numpy.linspace(-30, 100, 1000)
    alpha = 0.07
    ax.semilogy(phi, npfermi(phi, alpha), 'g')
    ax.semilogy(phi, dnpfermi(phi, alpha), 'g:')
    ax.semilogy(phi, num_npfermi(phi, alpha), 'g--', lw=2)

    phi = -5
    print fdk(0.5,phi)/_num_fermi(phi, alpha)
    print _npfermi(phi, alpha)/_num_npfermi(phi, alpha)
    phi = 5
    print fdk(0.5,phi)/_num_fermi(phi, alpha)
    print _npfermi(phi, alpha)/_num_npfermi(phi, alpha)
    phi = 30
    print fdk(0.5,phi)/_num_fermi(phi, alpha)
    print _npfermi(phi, alpha)/_num_npfermi(phi, alpha)
    plt.show()