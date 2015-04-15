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

from .units import units, to_units, cm, cm3
import numpy


# k = 8.6173324e-5 * unit.eV * unit.K**-1
# h = 4.135667516e-15 * unit.eV * unit.s
# pi = numpy.pi
# m0 = 9.10938215e-31 * unit.kg
# prefactor = (2*(2*pi*m0*k*unit.K / h**2)**(1.5)).to(1./cm3)
N_prefactor = 4.82936522463e15 # cm**-3 K**-1.5

class Material(object):
    def __init__(self, alloy, doping=None, Na=None, Nd=None):
        '''
        Arguments
        ---------
        alloy : openbandparams alloy
            Alloy.
        doping : float
            Ionized dopant concentration (cm**-3).
            Positive for p-type, negative for n-type.
        Na : float
            Ionized acceptor density (cm**-3)
        Nd : float
            Ionized donor density (cm**-3)
        '''
        # initialize caches
        self._Nc = {}
        self._Nv = {}
        self._nie = {}
        self._Ei = {}
        
        self._alloy = alloy

        if doping is not None and Na is not None:
            raise ValueError('If doping is specified, Na cannot be')
        if doping is not None and Nd is not None:
            raise ValueError('If doping is specified, Nd cannot be')

        if doping is not None:
            doping = to_units(doping, 1./cm3)
            if doping >= 0:
                self._Na = doping
                self._Nd = 0.
            else:
                self._Na = 0
                self._Nd = -doping
        else:
            if Na is not None:
                self._Na = to_units(Na, 1./cm3)
            else:
                self._Na = 0.
    
            if Nd is not None:
                self._Nd = to_units(Nd, 1./cm3)
            else:
                self._Nd = 0.
        # Nnet = Nd - Na
        self._Nnet = self._Nd - self._Na
        
            
    def __getattr__(self, name):
        if hasattr(self._alloy, name):
            # pass Material.Eg(), etc. on to the alloy
            return getattr(self._alloy, name)
        else:
            # Default behaviour
            raise AttributeError
    
    def Ei(self, T=300):
        if T in self._Ei:
            return self._Ei[T]
        else:
            k = 8.6173324e-5 # eV K**-1
            Ei = self.Eg()/2+k*T/2*numpy.log(self.Nv(T=T)/self.Nc(T=T))
            self._Ei[T] = Ei
            return Ei
    
    def nie(self, T=300):
        '''Returns the effective intrinsic carrier concentration (cm**-3)'''
        #TODO: include bandgap reduction, degeneracy, and non-parabolicity
        # using an additive energy that changes with doping.
        if T in self._nie:
            return self._nie[T]
        else:
            k = 8.6173324e-5 # eV K**-1
            Vt2 = k*T*2 # eV
            Eg = self.Eg(T=T)
            Nc = self.Nc(T=T)
            Nv = self.Nv(T=T)
            nie = numpy.sqrt(Nc*Nv)*numpy.exp(-Eg/Vt2)
            self._nie[T] = nie
            return nie
    
    def Nc(self, T=300):
        if T in self._Nc:
            return self._Nc[T]
        else:
            meff = self.meff_e_DOS(T=T)
            Nc = N_prefactor * (meff*T)**(1.5)
            self._Nc[T] = Nc
            return Nc

    def Nv(self, T=300):
        if T in self._Nv:
            return self._Nv[T]
        else:
            meff = .41#self.meff_h_DOS(T=T)
            Nv = N_prefactor * (meff*T)**(1.5)
            self._Nv[T] = Nv
            return Nv

    def meff_e_DOS(self, **kwargs):
        Eg_Gamma = self.Eg_Gamma(**kwargs)
        Eg_X = self.Eg_X(**kwargs)
        Eg_L = self.Eg_L(**kwargs)
        if Eg_Gamma < Eg_X:
            if Eg_Gamma < Eg_L:
                return self.meff_e_Gamma()
            else:
                return self.meff_e_L_DOS()
        else:
            if Eg_X < Eg_L:
                return self.meff_e_X_DOS()
            else:
                return self.meff_e_L_DOS()

    def meff_h_DOS(self, T=300):
        #TODO: calculate this properly
        return (self.meff_hh_100(T=T)*self.meff_hh_111(T=T))**0.5

    def Na(self, **kwargs):
        '''Returns the ionized acceptor concentration, Na (cm**-3)'''
        return self._Na

    def Nd(self, **kwargs):
        '''Returns the ionized donor concentration, Nd (cm**-3)'''
        return self._Nd

    def Nnet(self, **kwargs):
        '''Returns the net fixed charge, Nnet (cm**-3)'''
        return self._Nnet
