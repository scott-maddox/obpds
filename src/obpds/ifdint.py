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
'''
Fermi-Dirac integrals
    
[1] T. Fukushima, "Precise and fast computation of Fermi-Dirac integral
    of integer and half integer order by piecewise minimax rational
    approximation," Applied Mathematics and Computation, vol. 259,
    pp. 708-729, May 2015.
'''
import numpy
from numpy import exp, sqrt, log
from scipy.optimize import newton
from fdint import fdk, dfdk


__all__ = ['ifd1h']


def _ifd1h(nu):
    '''
    Inverse Fermi-Dirac integral of order 1/2.

    Parameters
    ----------
    nu : float
        normalized carrier concentration, n/Nc.
    
    Returns
    -------
    eta : float
        normalized Fermi energy, (phi_n-Ec)/Vt
    '''
    f = lambda eta: fdk(0.5,eta) - nu
    fprime = lambda eta: dfdk(0.5,eta)
    if nu < 10:
        guess = log(nu)
    else:
        guess = nu**1.5
    return newton(f, guess, fprime=fprime)
ifd1h = numpy.vectorize(_ifd1h)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    _, ax = plt.subplots()
    nu = numpy.logspace(-10, 3, 10000)
    ax.semilogx(nu, ifd1h(nu))
    plt.show()