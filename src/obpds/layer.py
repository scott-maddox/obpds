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

from .units import to_units, cm

class Layer(object):
    def __init__(self, thickness, material):
        '''
        Params
        ------
        thickness : float
            thickness (cm)
        material : openbandparams material
            the material in the layer
        '''
        self.thickness = to_units(thickness, cm)
        self.material = material

    def get_thickness(self):
        '''
        Returns the layer thickness (cm).
        '''
        return self.thickness

    def write_AMBER_recipe(self, fobj):
        fobj.write('l {} {:.3f} ! Angs\n'.format(self.material.name,
                                                 self.thickness*1e8))

    def get_flatband(self, T=300):
        '''
        returns x, Ev, Ec, Ei
        
        x will be [0, thickness]
        Ev will be [VBO, VBO]
        Ec will be [VBO+Eg, VBO+Eg]
        Ei will be [VBO+Ei, VBO+Ei]
        
        Arguments
        ---------
        T : float
            the temperature
        '''
        x = [0, self.get_thickness()]
        Ev = [self.material.VBO(T=T),
              self.material.VBO(T=T)]
        Ec = [self.material.VBO(T=T)+self.material.Eg(T=T),
              self.material.VBO(T=T)+self.material.Eg(T=T)]
        Ei = [self.material.VBO(T=T)+self.material.Ei(T=T),
              self.material.VBO(T=T)+self.material.Ei(T=T)]
        return x, Ev, Ec, Ei
        