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
from .solver import poisson_eq, poisson_zero_current
from .solution import ParametersSolution, FlatbandSolution


__all__ = ['TwoTerminalDevice']


class TwoTerminalDevice(object):
    '''
    A two terminal device composed of a number of layers with two contacts
    (left/top and right/bottom).
    '''
    def __init__(self, layers, contacts=None, Fp=None, Fn=None):
        '''
        Parameters
        ----------
        layers : list of `Layer`s
            layers
        contacts : list of `Contact`s (default=None)
            contacts; if None, defaults to two `OhmicContact`s
        Fp : str or list of str's
            Specifies control of the hole quasi-Fermi energy.
            'left' for control by the left contact.
            'right' for control by the right contact.
            None for Fp = inf.
        Fn : str or list of str's
            Specifies control of the electron quasi-Fermi energy.
            'left' for control by the left contact.
            'right' for control by the right contact.
            None for Fn = -inf.
        '''
        # Cache
        self._Fp_Fn = {}
        self._parameters = {}
        self._flatband = {}
        self._equilibrium = {}
        self._zero_current = {}
        
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
        
        if hasattr(Fp, '__iter__') and len(Fp) != len(layers):
            raise TypeError('len(Fp) != len(layers)')
        if hasattr(Fn, '__iter__') and len(Fn) != len(layers):
            raise TypeError('len(Fn) != len(layers)')
        self._Fp = Fp
        self._Fn = Fn
    
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
    
    def _calc_Fp_Fn(self, N):
        if self._Fp is None:
            Fp = [None]*N
        elif not hasattr(self._Fp, '__iter__'):
            Fp = [self._Fp]*N
        else:
            layer_xs = []
            last_x = 0
            for layer in self._layer:
                next_x = last_x + layer.get_thickness()
                layer_xs.append(next_x)
                last_x = next_x

            layer_xs = numpy.array(layer_xs)
            Fp = []
            for x in self._get_x(N):
                i = numpy.searchsorted(layer_xs, x)
                Fp.append(Fp[i])
                
        if self._Fn is None:
            Fn = [None]*N
        elif not hasattr(self._Fn, '__iter__'):
            Fn = [self._Fn]*N
        else:
            layer_xs = []
            last_x = 0
            for layer in self._layer:
                next_x = last_x + layer.get_thickness()
                layer_xs.append(next_x)
                last_x = next_x

            layer_xs = numpy.array(layer_xs)
            Fn = []
            for x in self._get_x(N):
                i = numpy.searchsorted(layer_xs, x)
                Fn.append(Fn[i])
        
        return Fp, Fn

    def _get_Fp_Fn(self, N):
        if N in self._Fp_Fn:
            return self._Fp_Fn[N]
        else:
            s = self._calc_Fp_Fn(N)
            self._Fp_Fn[N] = s
            return s

    def _get_materials(self, N):
        return [self._layer.get_material(x_i) for x_i in self._get_x(N)]

    def _calc_parameters(self, T, N):
        s = ParametersSolution()
        s.materials = materials = self._get_materials(N)
        for p in ['VBO', 'CBO_Gamma', 'CBO_L', 'CBO_X', 'CBO', 'Ei', 'ni',
                  'dielectric', 'Na', 'Nd', 'Nnet', 'Nc_Gamma', 'Nc_L',
                  'Nc_X', 'Nc', 'Nv', 'nonparabolicity', 'electron_affinity']:
            value = numpy.array([getattr(m, p)(T=T) for m in materials],
                                dtype=float)
            setattr(s, p, value)
        return s
    
    def _get_parameters(self, T, N):
        if (T, N) in self._parameters:
            return self._parameters[(T, N)]
        else:
            s = self._calc_parameters(T, N)
            self._parameters[(T, N)] = s
            return s

    def _calc_flatband(self, T, N):
        x = self._get_x(N)
        p = self._get_parameters(T, N)
        return FlatbandSolution(T=T, N=N, x=x,
                                Ev=p.VBO,
                                Ec_Gamma=p.CBO_Gamma,
                                Ec_L=p.CBO_L,
                                Ec_X=p.CBO_X,
                                Ec=p.CBO,
                                Ei=p.VBO+p.Ei)
    
    def _get_flatband(self, T, N):
        if (T, N) in self._flatband:
            return self._flatband[(T, N)]
        else:
            s = self._calc_flatband(T, N)
            self._flatband[(T, N)] = s
            return s
    
    def _calc_equilibrium(self, T, N, approx='parabolic'):
        solution = poisson_eq(self, T=T, N=N, approx=approx)
        self._equilibrium[(T, N, approx)] = solution
        return solution
    
    def get_equilibrium(self, T=300., N=1000, approx='parabolic'):
        '''
        Returns an `EquilibriumSolution` instance.
        
        Arguments
        ---------
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default='parabolic')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        if (T, N, approx) in self._equilibrium:
            return self._equilibrium[(T, N, approx)]
        else:
            return self._calc_equilibrium(T, N, approx)

    def show_equilibrium(self, T=300., N=1000, approx='parabolic'):
        '''
        Plot and show the band profile at equilibrium.
        
        Arguments
        ---------
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default='parabolic')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        solution = self.get_equilibrium(T, N, approx)
        x = solution.x*1e7 # nm
        import matplotlib.pyplot as plt
        _, (ax1, ax2) = plt.subplots(2, 1, sharex='col')
        ax1.set_ymargin(0.05)
        ax2.set_ymargin(0.05)
        ax1.plot(x, solution.Ev, 'r-', label='$E_v$')
        ax1.plot(x, solution.Ec, 'b-', label='$E_c$')
        ax1.plot(x, solution.Ef, 'k--', label='$E_f$')
        ax1.plot(x, solution.Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        if (solution.Na > 0.).any():
            ax2.semilogy(x, solution.Na, 'r-', label='$N_A$')
        if (solution.Nd > 0.).any():
            ax2.semilogy(x, solution.Nd, 'b-', label='$N_D$')
        ax2.semilogy(x, solution.p, 'r--', label='$p$')
        ax2.semilogy(x, solution.n, 'b--', label='$n$')
        ax2.set_ylabel('Concentration (cm$^{-3}$)')
        ax2.set_xlabel('Depth (nm)')
        ymin, ymax = ax2.get_ylim()
        if ymax/ymin > 1e10:
            ax2.set_ylim(ymax/1e10, ymax)
        plt.show()

    def save_equilibrium(self, path, show=False, T=300, N=1000, approx='parabolic'):
        '''
        Save the bands at equilibrium.
        
        Arguments
        ---------
        path : string
            the file path
        show : bool
            shows the bands if True
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default='parabolic')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        if show:
            self.show_equilibrium(T, N, approx)
        s = self.get_equilibrium(T, N, approx)
        with open(path, 'w') as f:
            f.write('x\tEv\tEc\tEi\tp\tn\tNa\tNd\n')
            for i in xrange(s.x.size):
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                        ''.format(s.x[i], s.Ev[i], s.Ec[i], s.Ei[i],
                                  s.p[i], s.n[i], s.Na[i], s.Nd[i]))
    
    def _calc_zero_current(self, V, T, N, approx='parabolic'):
        solution = poisson_zero_current(self, V=V, T=T, N=N, approx=approx)
        self._zero_current[(V, T, N, approx)] = solution
        return solution
    
    def get_zero_current(self, V, T=300., N=1000, approx='parabolic'):
        '''
        Returns a `ZeroCurrentSolution` instance.
        
        Arguments
        ---------
        V : float
            Bias voltage, i.e. left/top contact bias - right/bottom contact bias
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default='parabolic')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        if (V, T, N, approx) in self._zero_current:
            return self._zero_current[(V, T, N, approx)]
        else:
            return self._calc_zero_current(V, T, N, approx)

    def show_zero_current(self, V, T=300., N=1000, approx='parabolic'):
        '''
        Plot and show the band profile at a given bias voltage under the
        zero-current approximation.
        
        Arguments
        ---------
        V : float
            Bias voltage, i.e. left/top contact bias - right/bottom contact bias
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default='parabolic')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        solution = self.get_zero_current(V, T, N, approx)
        x = solution.x*1e7 # nm
        import matplotlib.pyplot as plt
        _, (ax1, ax2) = plt.subplots(2, 1, sharex='col')
        ax1.set_ymargin(0.05)
        ax2.set_ymargin(0.05)
        ax1.plot(x, solution.Ev, 'r-', label='$E_v$')
        ax1.plot(x, solution.Ec, 'b-', label='$E_c$')
        ax1.plot(x, solution.Fp, 'r--', label='$F_p$')
        ax1.plot(x, solution.Fn, 'b--', label='$F_n$')
        ax1.plot(x, solution.Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        if (solution.Na > 0.).any():
            ax2.semilogy(x, solution.Na, 'r-', label='$N_A$')
        if (solution.Nd > 0.).any():
            ax2.semilogy(x, solution.Nd, 'b-', label='$N_D$')
        ax2.semilogy(x, solution.p, 'r--', label='$p$')
        ax2.semilogy(x, solution.n, 'b--', label='$n$')
        ax2.set_ylabel('Concentration (cm$^{-3}$)')
        ax2.set_xlabel('Depth (nm)')
        ymin, ymax = ax2.get_ylim()
        if ymax/ymin > 1e10:
            ax2.set_ylim(ymax/1e10, ymax)
        plt.show()