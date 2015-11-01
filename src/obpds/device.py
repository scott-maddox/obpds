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
from .solver import (poisson_eq, poisson_zero_current,
                     capacitance_zero_current)
from .solution import ParametersSolution, FlatbandSolution
from .config import cfg


__all__ = ['TwoTerminalDevice']

BLUE = '#072a9f'
RED = '#d62315'

# electron charge
q = 1.602176565e-19 # C

#TODO: test and use this decorator
def _cached_method(f):
    cache_name = '__%s_cache' % f.__name__
    def wrapper(self, *args):
        cache = getattr(self, cache_name, None)
        if cache is None:
            cache = {}
            setattr(self, cache_name, cache)
        if args in cache:
            return cache[args]
        res = cache[args] = f(self, *args)
        return res
    return wrapper

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
        self._capacitance = {}
        
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

    #TODO: make this consistent with the other show_* methods.
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
                Fp.append(self._Fp[i])
                
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
                Fn.append(self._Fn[i])
        
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
    
    def _calc_equilibrium(self, T, N, approx='kane'):
        solution = poisson_eq(self, T=T, N=N, approx=approx)
        self._equilibrium[(T, N, approx)] = solution
        return solution
    
    def get_thickness(self):
        return self._layer.get_thickness()
    
    def _has_equilibrium(self, T=300., N=1000, approx='kane'):
        '''
        Returns True if the equilbrium solution is cached.
        
        Arguments
        ---------
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        return (T, N, approx) in self._equilibrium

    def get_equilibrium(self, T=300., N=1000, approx='kane'):
        '''
        Returns an `EquilibriumSolution` instance.
        
        Arguments
        ---------
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        if self._has_equilibrium(T, N, approx):
            return self._equilibrium[(T, N, approx)]
        else:
            return self._calc_equilibrium(T, N, approx)

    def show_equilibrium(self, T=300., N=1000, approx='kane'):
        '''
        Plot and show the band profile at equilibrium.
        
        Arguments
        ---------
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        solution = self.get_equilibrium(T, N, approx)
        x = solution.x*1e7 # nm
        import matplotlib.pyplot as plt
        plt.style.use(['ggplot'])
        _, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col',
                                          figsize=(10, 10),
                                          tight_layout=cfg['plot/tight_layout'])
        ax1.set_ymargin(0.05)
        ax1.plot(x, solution.Ev, 'r-', label='$E_v$')
        ax1.plot(x, solution.Ec, 'b-', label='$E_c$')
        ax1.plot(x, solution.Ef, 'k--', label='$E_f$')
        ax1.plot(x, solution.Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        
        ax2.set_ymargin(0.05)
        if (solution.Na > 0.).any():
            ax2.semilogy(x, solution.Na, 'r-', label='$N_A$')
        if (solution.Nd > 0.).any():
            ax2.semilogy(x, solution.Nd, 'b-', label='$N_D$')
        ax2.semilogy(x, solution.p, 'r--', label='$p$')
        ax2.semilogy(x, solution.n, 'b--', label='$n$')
        ax2.set_ylabel('Concentration (cm$^{-3}$)')
        ymin, ymax = ax2.get_ylim()
        if ymax/ymin > cfg['plot/semilogy/yrange']:
            ax2.set_ylim(ymax/cfg['plot/semilogy/yrange'], ymax)
            
        ax3.set_ymargin(0.05)
        (ax3_field,) = ax3.plot(x, solution.field, 'k-')
        (ax3_dEv_dx,) = ax3.plot(x, solution.dEv_dx, 'r-', alpha=0.5)
        (ax3_dEc_dx,) = ax3.plot(x, solution.dEc_dx, 'b-', alpha=0.5)
        ax3.set_ylabel('Effective Field (V/cm)')
        ax3.set_xlabel('Depth (nm)')
        ax3.yaxis.get_major_formatter().set_powerlimits((-3, 3))
        self.filtered_autolim(ax3, solution.dEv_dx, solution.dEc_dx,
                              solution.field)
        
        plt.show()
    
    def _save_solution(self, s, path):
        names = [name for name in dir(s) if not name.startswith('_')]
        excludes = ['V', 'N', 'T', 'materials']
        for exclude in excludes:
            if exclude in names:
                names.remove(exclude)
        arrays = [getattr(s, name) for name in names]
        if not names:
            return
        header = '\t'.join(names)+'\n'
        template = '\t'.join(['{}' for name in names])+'\n'
        with open(path, 'w') as f:
            f.write(header)
            for i in xrange(arrays[0].size):
                values = [repr(array[i]) for array in arrays]
                f.write(template.format(*values))

    def save_equilibrium(self, path, show=False, T=300, N=1000, approx='kane'):
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
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        if show:
            self.show_equilibrium(T, N, approx)
        s = self.get_equilibrium(T, N, approx)
        self._save_solution(s, path)
    
    def _calc_zero_current(self, V, T, N, approx):
        return poisson_zero_current(self, V=V, T=T, N=N, approx=approx)
    
    def has_zero_current(self, V, T=300., N=1000, approx='kane'):
        '''
        Returns True if the zero current solution is cached.
        
        Arguments
        ---------
        V : float
            Bias voltage, i.e. left/top contact bias - right/bottom contact bias
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        return (V, T, N, approx) in self._zero_current
    
    def get_zero_current(self, V, T=300., N=1000, approx='kane'):
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
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        if self.has_zero_current(V, T, N, approx):
            return self._zero_current[(V, T, N, approx)]
        else:
            solution = self._calc_zero_current(V, T, N, approx)
            self._zero_current[(V, T, N, approx)] = solution
            return solution

    def show_zero_current(self, V, T=300., N=1000, approx='kane'):
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
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        self._zero_current_image(V, path=None, show=True,
                                 T=300., N=1000, approx='kane')

    def save_zero_current_image(self, path, V, show=False,
                                T=300., N=1000, approx='kane'):
        self._zero_current_image(V=V, path=path, show=show,
                                 T=300., N=1000, approx='kane')

    def _zero_current_image(self, V, path=None, show=False, 
                            T=300., N=1000, approx='kane'):
        solution = self.get_zero_current(V, T, N, approx)
        x = solution.x*1e7 # nm
        import matplotlib.pyplot as plt
        plt.style.use(['ggplot'])
        _, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col',
                                          figsize=(10, 10),
                                          tight_layout=cfg['plot/tight_layout'])

        ax1.set_ymargin(0.05)
        ax2.set_ymargin(0.05)
        ax1.plot(x, solution.Ev, 'r-', label='$E_v$')
        ax1.plot(x, solution.Ec, 'b-', label='$E_c$')
        ax1.plot(x, solution.Fp, 'r--', label='$F_p$')
        ax1.plot(x, solution.Fn, 'b--', label='$F_n$')
        ax1.plot(x, solution.Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        
        ax2.set_ymargin(0.05)
        if (solution.Na > 0.).any():
            ax2.semilogy(x, solution.Na, 'r-', label='$N_A$')
        if (solution.Nd > 0.).any():
            ax2.semilogy(x, solution.Nd, 'b-', label='$N_D$')
        ax2.semilogy(x, solution.p, 'r--', label='$p$')
        ax2.semilogy(x, solution.n, 'b--', label='$n$')
        ax2.set_ylabel('Concentration (cm$^{-3}$)')
        ymin, ymax = ax2.get_ylim()
        if ymax/ymin > cfg['plot/semilogy/yrange']:
            ax2.set_ylim(ymax/cfg['plot/semilogy/yrange'], ymax)
            
        ax3.set_ymargin(0.05)
        dEv_dx = numpy.empty_like(solution.dEv_dx)
        dEc_dx = numpy.empty_like(solution.dEc_dx)
        (ax3_field,) = ax3.plot(x, solution.field, 'k-')
        (ax3_dEv_dx,) = ax3.plot(x, solution.dEv_dx, 'r-', alpha=0.5)
        (ax3_dEc_dx,) = ax3.plot(x, solution.dEc_dx, 'b-', alpha=0.5)
        ax3.axhline(0, color='grey')
        ax3.set_ylabel('Effective Field (V/cm)')
        ax3.set_xlabel('Depth (nm)')
        ax3.yaxis.get_major_formatter().set_powerlimits((-3, 3))
        self.filtered_autolim(ax3, solution.dEv_dx, solution.dEc_dx,
                              solution.field)
        if path is not None:
            plt.savefig(path)
        if show:
            plt.show()

    def interactive_zero_current(self, T=300., N=1000, approx='kane'):
        '''
        Arguments
        ---------
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        solution = self.get_zero_current(0., T, N, approx)
        x = solution.x*1e7 # nm
        import matplotlib.pyplot as plt
        plt.style.use(['ggplot'])
        fig = plt.figure(figsize=(10, 10),
                         #facecolor='white', edgecolor='white',
                         tight_layout=cfg['plot/tight_layout'])
        ax3 = plt.subplot2grid(shape=(3,15), loc=(2,0), colspan=14)
        ax2 = plt.subplot2grid(shape=(3,15), loc=(1,0), colspan=14, sharex=ax3)
        ax1 = plt.subplot2grid(shape=(3,15), loc=(0,0), colspan=14, sharex=ax3)
        ax4 = plt.subplot2grid(shape=(3,15), loc=(0,14), sharey=ax1)
        ax4.get_xaxis().set_visible(False)
        ax4.get_yaxis().set_visible(False)
        ax4.axhline(0, color='grey')
        ax4_V = ax4.axhline(y=0., color='k', linestyle='-')

        ax1.set_ymargin(0.1)
        ax2.set_ymargin(0.1)
        ax3.set_ymargin(0.1)

        (ax1_Ev,) = ax1.plot(x, solution.Ev, 'r-', label='$E_v$')
        (ax1_Ec,) = ax1.plot(x, solution.Ec, 'b-', label='$E_c$')
        (ax1_Fp,) = ax1.plot(x, solution.Fp, 'r--', label='$F_p$')
        (ax1_Fn,) = ax1.plot(x, solution.Fn, 'b--', label='$F_n$')
        (ax1_Ei,) = ax1.plot(x, solution.Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        
        nans = numpy.empty_like(x)
        nans.fill(numpy.nan)
        (ax2_Na,) = ax2.plot(x, solution.Na, 'r-', label='$N_A$')
        (ax2_Nd,) = ax2.plot(x, solution.Nd, 'b-', label='$N_D$')
        (ax2_p,) = ax2.semilogy(x, solution.p, 'r--', label='$p$')
        (ax2_n,) = ax2.semilogy(x, solution.n, 'b--', label='$n$')
        ax2.set_ylabel('Concentration (cm$^{-3}$)')
        ymin, ymax = ax2.get_ylim()
        if ymax/ymin > cfg['plot/semilogy/yrange']:
            ax2.set_ylim(ymax/cfg['plot/semilogy/yrange'], ymax)
        
        (ax3_field,) = ax3.plot(x, solution.field, 'k-')
        (ax3_dEv_dx,) = ax3.plot(x, solution.dEv_dx, color=RED, alpha=0.7)
        (ax3_dEc_dx,) = ax3.plot(x, solution.dEc_dx, color=BLUE, alpha=0.7)
        ax3.axhline(0, color='grey')
        ax3.set_ylabel('Effective Field (V/cm)')
        ax3.set_xlabel('Depth (nm)')
        ax3.yaxis.get_major_formatter().set_powerlimits((-3, 3))
        self.filtered_autolim(ax3, solution.dEv_dx, solution.dEc_dx,
                              solution.field)

#         new_ymin = solution.field.min()
#         new_ymax = solution.field.max()
#         if ax3._ymargin > 0:
#             delta = (new_ymax - new_ymin) * ax3._ymargin
#             new_ymin -= delta
#             new_ymax += delta
#         ax3.set_ybound(new_ymin, new_ymax)
        
        def onclick(event):
            if event.inaxes != ax4:
                return
            V = event.ydata
            ax4_V.set_data(([0, 1], [V, V]))
            solution = self.get_zero_current(V, T, N, approx)
            ax1_Ev.set_data(x, solution.Ev)
            ax1_Ec.set_data(x, solution.Ec)
            ax1_Fp.set_data(x, solution.Fp)
            ax1_Fn.set_data(x, solution.Fn)
            ax1_Ei.set_data(x, solution.Ei)
#             if (solution.Na > 0.).any():
#                 ax2_Na.set_data(x, solution.Na)
#             else:
#                 ax2_Na.set_data(x, nans)
#             if (solution.Nd > 0.).any():
#                 ax2_Nd.set_data(x, solution.Nd)
#             else:
#                 ax2_Nd.set_data(x, solution.Nd)
            ax2_p.set_data(x, solution.p)
            ax2_n.set_data(x, solution.n)
            ax3_dEv_dx.set_data(x, solution.dEv_dx)
            ax3_dEc_dx.set_data(x, solution.dEc_dx)
            ax3_field.set_data(x, solution.field)
            
            for ax in [ax1, ax2]:
                old_ymin, old_ymax = ax.get_ylim()
                ax.relim()
                new_ymin, new_ymax = ax.dataLim.intervaly
                if ax._ymargin > 0:
                    delta = (new_ymax - new_ymin) * ax._ymargin
                    new_ymin -= delta
                    new_ymax += delta
                ax.set_ybound(min(old_ymin, new_ymin), max(old_ymax, new_ymax))
                
            ymin, ymax = ax2.get_ylim()
            if ymax/ymin > cfg['plot/semilogy/yrange']:
                ax2.set_ybound(ymax/cfg['plot/semilogy/yrange'], ymax)
            
            self.filtered_autolim(ax3, solution.dEv_dx, solution.dEc_dx,
                                  solution.field)
            
            fig.canvas.draw()

        _cid = fig.canvas.mpl_connect('button_press_event', onclick)
        
        plt.show()
    
    @classmethod
    def filtered_autolim(cls, ax, *fields):
        threshold1 = 10
        threshold2 = 2
        fmin = numpy.inf
        fmax = -numpy.inf
        for field in fields:
            for i in xrange(2, field.size-2):
                delta_m2 = abs(field[i]-field[i-2])
                delta_m1 = abs(field[i]-field[i-1])
                delta_p1 = abs(field[i]-field[i+1])
                delta_p2 = abs(field[i]-field[i+2])
    #             if (delta_m1 ~= delta_m2 and
    #                 delta_p1 << delta_m1 and
    #                 delta_p2 >> delta_p1):
                if (delta_m1            < delta_m2*threshold2   and
                    delta_p1*threshold1 < delta_m1              and
                    delta_p2            > delta_p1*threshold1):
                    continue
    #             if (delta_m1 << delta_m2 and
    #                 delta_p1 >> delta_m1 and
    #                 delta_p2 ~= delta_p1):
                if (delta_m1*threshold1 < delta_m2            and
                    delta_p1            > delta_m1*threshold1 and
                    delta_p2            < delta_p1*threshold2):
                    continue
                fmin = min(fmin, field[i])
                fmax = max(fmax, field[i])
        if ax._ymargin > 0:
            delta = (fmax - fmin) * ax._ymargin
            fmin -= delta
            fmax += delta
        ax.set_ybound(fmin, fmax)

    def save_zero_current(self, V, path, show=False, T=300, N=1000,
                          approx='kane'):
        '''
        Save the band profile data and image at a given bias voltage under the
        zero-current approximation.
        
        Arguments
        ---------
        V : float
            Bias voltage, i.e. left/top contact bias - right/bottom contact bias
        path : string
            the file path without (excluding file extension)
        show : bool
            shows the bands if True
        T : float (default=300.)
            Device temperature
        N : int (default=1000)
            Number of grid points
        approx : str (default ='kane')
            If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
            bands approximation (fastest). If 'parabolic', use the parabolic
            bands approximation (fast). If 'kane', include Gamma-valley
            non-parabolicity under the k.p Kane approximation (slow).
        '''
        s = self.get_zero_current(V, T, N, approx)
        self._save_solution(s, path=path+'.txt')
        self._zero_current_image(V=V, path=path+'.png', show=show,
                                 T=T, N=N, approx=approx)
    
    def _calc_capacitance(self, V, dV, T=300, N=1000, approx='kane'):
        return capacitance_zero_current(self, V, dV, T, N, approx)
    
    def get_capacitance(self, V, dV=1e-3, T=300, N=1000, approx='kane'):
        '''
        Returns
        -------
        C : float
            capacitance in units of F/cm**2
        '''
        if (V, dV, T, N, approx) in self._capacitance:
            return self._capacitance[(V, dV, T, N, approx)]
        else:
            s = self._calc_capacitance(V, dV, T, N, approx)
            self._capacitance[(V, dV, T, N, approx)] = s
            return s
    
    def get_cv(self, Vstart, Vstop, Vnum=100, dV=1e-3, T=300, N=1000,
               approx='kane'):
        '''
        Returns
        -------
        C : ndarray
            capacitance in units of F/cm**2
        V : ndarray
            bias voltage in units of V
        '''
        V = numpy.linspace(Vstart, Vstop, Vnum)
        C = numpy.empty(Vnum)
        for i in xrange(Vnum):
            C[i] = self.get_capacitance(V[i], dV, T, N, approx)
        return C, V

    def show_cv(self, Vstart, Vstop, Vnum=50, dV=1e-3, T=300, N=1000,
                approx='kane'):
        C, V = self.get_cv(Vstart, Vstop, Vnum, dV, T, N, approx)
        rCs = 1/C**2
        ndV_drCs = (-1.)/numpy.gradient(rCs, (V[1]-V[0]))
        import matplotlib.pyplot as plt
        plt.style.use(['ggplot'])
        _, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col',
                                          figsize=(10, 10),
                                          tight_layout=cfg['plot/tight_layout'])
            
        ax1.set_ymargin(0.05)
        ax1.plot(V, C, 'r-')
#         ax1.axhline(0, color='grey')
        ax1.set_ylabel('Capacitance (F/cm$^2$)')
        
        ax2.set_ymargin(0.05)
        ax2.plot(V, rCs, 'r-')
        ax2.set_ylabel('1/C$^2$ (cm$^4$/F$^2$)')
        
        ax3.set_ymargin(0.05)
        try:
            ax3.semilogy(V, ndV_drCs, 'r-')
        except:
            ax3.set_yscale('linear')
#             ax3.plot(V, dV/numpy.gradient(1/C**2), 'r-')
        ax3.set_ylabel('-dV/d(1/C$^2$)')
        ax3.set_xlabel('Bias (V)')
        
        plt.show()
    
    def save_cv(self, path, Vstart, Vstop, Vnum=50, dV=1e-3, show=False,
                T=300, N=1000, approx='kane'):
        C, V = self.get_cv(Vstart, Vstop, Vnum, dV, T, N, approx)
        rCs = 1/C**2
        ndV_drCs = (-1.)/numpy.gradient(rCs, (V[1]-V[0]))
        if show:
            import matplotlib.pyplot as plt
            plt.style.use(['ggplot'])
            _, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col',
                                              figsize=(10, 10),
                                              tight_layout=cfg['plot/tight_layout'])
                
            ax1.set_ymargin(0.05)
            ax1.plot(V, C, 'r-')
    #         ax1.axhline(0, color='grey')
            ax1.set_ylabel('Capacitance (F/cm$^2$)')
            
            ax2.set_ymargin(0.05)
            ax2.plot(V, rCs, 'r-')
            ax2.set_ylabel('1/C$^2$ (cm$^4$/F$^2$)')
            
            ax3.set_ymargin(0.05)
            try:
                ax3.semilogy(V, ndV_drCs, 'r-')
            except:
                ax3.set_yscale('linear')
    #             ax3.plot(V, dV/numpy.gradient(1/C**2), 'r-')
            ax3.set_ylabel('-dV/d(1/C$^2$)')
            ax3.set_xlabel('Bias (V)')
            
            plt.show()
        
        header = 'V\tC/A\t1/C^2\tdV/d(1/C^2)\n'
        template = '{V}\t{C}\t{rCs}\t{ndV_drCs}\n'
        with open(path, 'w') as f:
            f.write(header)
            for i in xrange(V.size):
                f.write(template.format(V=V[i], C=C[i],
                                        rCs=rCs[i], ndV_drCs=ndV_drCs[i]))