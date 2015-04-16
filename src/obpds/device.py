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

from .layer import CompoundLayer
from .contact import Contact, OhmicContact
from .solver import poisson_eq
from .solution import FlatbandSolution


__all__ = ['TwoTerminalDevice']


class TwoTerminalDevice(object):
    '''
    A two terminal device composed of a number of layers with two contacts
    (left/top and right/bottom).
    '''
    def __init__(self, layers, contacts=None):
        '''
        Parameters
        ----------
        layers : list of `Layer`s
            layers
        contacts : list of `Contact`s (default=None)
            contacts; if None, defaults to two `OhmicContact`s
        '''
        # Cache
        self._equilibrium = {}
        self._flatband = {}
        
        self._layer = CompoundLayer(layers)

        if contacts is None:
            self._contacts = [OhmicContact(), OhmicContact()]
        elif len(contacts) != 2:
            raise ValueError('There must be exactly two contacts.')
        else:
            for contact in contacts:
                if not isinstance(contact, Contact):
                    raise TypeError('Contacts must be instances of '
                                    'the `Contact` class.')
            self._contacts = contacts
    
    def get_flatband(self, T=300.):
        '''
        returns x, Ev, Ec, Ei
        
        x will be numpy.array([0, ..., thickness])
        Ev will be numpy.array([VBO, ..., VBO])
        Ec will be numpy.array([CBO, ..., CBO])
        Ei will be numpy.array([VBO+Ei, ..., VBO+Ei])
        
        Arguments
        ---------
        T : float
            the temperature
        '''
        x, Ev, Ec, Ei = self._layer.get_flatband(T)
        return numpy.array(x), numpy.array(Ev), numpy.array(Ec), numpy.array(Ei)

    def show_flatband(self, T=300.):
        '''
        Show a plot of the band profile at flatband.
        
        Arguments
        ---------
        T : float (default=300.)
            the temperature
        '''
        import matplotlib.pyplot as plt
        _, ax = plt.subplots()
        x, Ev, Ec, Ei = self.get_flatband(T=T)
        x = x*1e7 # nm
        ax.plot(x, Ev, 'r-', label='$E_v$')
        ax.plot(x, Ec, 'b-', label='$E_c$')
        ax.plot(x, Ei, 'k:', label='$E_i$')
        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel('Depth (nm)')
        plt.show()
    
    def _get_x(self, N):
        return numpy.linspace(0, self._layer.get_thickness(), N)

    def _get_materials(self, N):
        return [self._layer.get_material(x_i) for x_i in self._get_x(N)]

    def _calc_flatband(self, T, N):
        x = self._get_x(N)
        materials = self._get_materials(N)
        Ev = numpy.array([m.VBO() for m in materials], dtype=float)
        Ec = numpy.array([m.CBO() for m in materials], dtype=float)
        Ei = numpy.array([m.VBO()+m.Ei() for m in materials], dtype=float)
        solution = FlatbandSolution(T, N, x, Ev, Ec, Ei)
        self._flatband[(T, N)] = solution
        return solution
    
    def _get_flatband(self, T, N):
        if (T, N) in self._flatband:
            return self._flatband[(T, N)]
        else:
            return self._calc_flatband(T, N)
    
    def _calc_equilibrium(self, T, N):
        solution = poisson_eq(self, T=T, N=N)
        self._equilibrium[(T, N)] = solution
        return solution
    
    def get_equilibrium(self, T=300., N=1000):
        '''
        Returns an `EquilibriumSolution` instance.
        '''
        if (T, N) in self._equilibrium:
            return self._equilibrium[(T, N)]
        else:
            return self._calc_equilibrium(T, N)

    def show_equilibrium(self, T=300., N=1000):
        '''
        Show a plot of the band profile at equilibrium.
        
        Arguments
        ---------
        T : float
            the temperature
        N : int
            the number of grid points
        '''
        solution = self.get_equilibrium(T, N)
        x = solution.x*1e7 # nm
        import matplotlib.pyplot as plt
        _, (ax1, ax2) = plt.subplots(2, 1, sharex='col')
        ax1.plot(x, solution.Ev, 'r-', label='$E_v$')
        ax1.plot(x, solution.Ec, 'b-', label='$E_c$')
        ax1.plot(x, solution.Ef, 'k--', label='$E_f$')
        ax1.plot(x, solution.Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        ax2.semilogy(x, solution.Na, 'r-', label='$N_A$')
        ax2.semilogy(x, solution.Nd, 'b-', label='$N_D$')
        ax2.semilogy(x, solution.p, 'r--', label='$p$')
        ax2.semilogy(x, solution.n, 'b--', label='$n$')
        ax2.set_ylabel('Concentration (cm$^{-3}$)')
        ax2.set_xlabel('Depth (nm)')
        plt.show()

    def save_equilibrium(self, path, show=False, T=300, N=1000):
        '''
        Save the bands at equilibrium.
        
        Arguments
        ---------
        path : string
            the file path
        show : bool
            shows the bands if True
        T : float
            the temperature
        N : int
            the number of grid points
        '''
        if show:
            self.show_equilibrium(T=T, N=N)
        x, Ev, Ec, Ei, p, n, Na, Nd = poisson_eq(self, T=T, N=N)
        with open(path, 'w') as f:
            f.write('x\tEv\tEc\tEi\tp\tn\tNa\tNd\n')
            for i in xrange(x.size):
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                        ''.format(x[i], Ev[i], Ec[i], Ei[i],
                                  p[i], n[i], Na[i], Nd[i]))