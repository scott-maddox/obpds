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


__all__ = ['Material']


# Boltzmann constant
k = 8.6173324e-5 # eV K**-1
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
        self._Nc_Gamma = {}
        self._Nc_X = {}
        self._Nc_L = {}
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
            Ei = self.Eg()/2.+k*T/2.*numpy.log(self.Nv(T=T)/self.Nc(T=T))
            self._Ei[T] = Ei
            return Ei
    
    def ni(self, T=300):
        return (numpy.sqrt(self.Nv(T=T)*self.Nc(T=T)) *
                numpy.exp(-self.Eg(T=T)/(2*k*T)))
    
    def Nc_Gamma(self, T=300):
        if T in self._Nc_Gamma:
            return self._Nc_Gamma[T]
        else:
            meff = self.meff_e_Gamma(T=T)
            Nc = N_prefactor * (meff*T)**(1.5)
            self._Nc_Gamma[T] = Nc
            return Nc
    
    def Nc_L(self, T=300):
        if T in self._Nc_L:
            return self._Nc_L[T]
        else:
            meff = self.meff_e_L_DOS(T=T)
            Nc = N_prefactor * (meff*T)**(1.5)
            self._Nc_L[T] = Nc
            return Nc
    
    def Nc_X(self, T=300):
        if T in self._Nc_X:
            return self._Nc_X[T]
        else:
            meff = self.meff_e_X_DOS(T=T)
            Nc = N_prefactor * (meff*T)**(1.5)
            self._Nc_X[T] = Nc
            return Nc
    
    def Nc(self, T=300.):
        Eg_Gamma = self.Eg_Gamma(T=T)
        Eg_X = self.Eg_X(T=T)
        Eg_L = self.Eg_L(T=T)
        if Eg_Gamma < Eg_X:
            if Eg_Gamma < Eg_L:
                return self.Nc_Gamma()
            else:
                return self.Nc_L()
        else:
            if Eg_X < Eg_L:
                return self.Nc_X()
            else:
                return self.Nc_L()

    def Nv(self, T=300.):
        if T in self._Nv:
            return self._Nv[T]
        else:
            meff = self.meff_h_DOS(T=T)
            Nv = N_prefactor * (meff*T)**(1.5)
            self._Nv[T] = Nv
            return Nv

    def meff_h_DOS(self, T=300.):
        #TODO: calculate this properly
        return self.meff_hh_100(T=T)

    def Na(self, **kwargs):
        '''Returns the ionized acceptor concentration, Na (cm**-3)'''
        return self._Na

    def Nd(self, **kwargs):
        '''Returns the ionized donor concentration, Nd (cm**-3)'''
        return self._Nd

    def Nnet(self, **kwargs):
        '''Returns the net fixed charge, Nnet (cm**-3)'''
        return self._Nnet
