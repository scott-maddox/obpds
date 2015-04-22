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
        material : Material
            material the layer is composed of
        '''
        self._thickness = to_units(thickness, cm)
        self._material = material

    def get_thickness(self):
        '''
        Returns the layer thickness (cm).
        '''
        return self._thickness

    def get_material(self, x):
        '''
        Returns the material at a given position (cm) within the layer.
        
        Parameters
        ----------
        x : float
            position (cm) between zero and the layer thickness
        '''
        if x < 0. or x > self.get_thickness():
            raise ValueError('x not within range [{:g}, {:g}]'
                             ''.format(0., self.get_thickness()))
        return self._material

    def get_flatband(self, T=300.):
        '''
        returns x, Ev, Ec, Ei
        
        x will be [0, thickness]
        Ev will be [VBO, VBO]
        Ec will be [CBO, CBO]
        Ei will be [VBO+Ei, VBO+Ei]
        
        Arguments
        ---------
        T : float
            the temperature
        '''
        x = (0., self.get_thickness())
        VBO = self._material.VBO(T=T)
        CBO = self._material.CBO(T=T)
        ILO = VBO + self._material.Ei(T=T)
        return x, (VBO, VBO), (CBO, CBO), (ILO, ILO)

    def write_AMBER_recipe(self, fobj):
        fobj.write('l {} {:.3f} ! Angs\n'.format(self._material.name,
                                                 self._thickness*1e8))

class CompoundLayer(Layer):
    
    def __init__(self, layers):
        self._layers = []
        if len(layers) < 1:
            raise ValueError('There must be at least one layer.')
        else:
            for layer in layers:
                self.append(layer)

    def __iter__(self):
        return self._layers.__iter__()
    
    def append(self, layer):
        '''
        Append a layer on the right/bottom.
        '''
        if not isinstance(layer, Layer):
            raise TypeError('Layers must be an instance of '
                            'the `Layer` class.')
        self._layers.append(layer)
    
    def insert(self, index, layer):
        '''
        Insert a layer at the given index.
        '''
        if not isinstance(layer, Layer):
            raise TypeError('Layers must be an instance of '
                            'the `Layer` class.')
        self._layers.insert(index, layer)

    def get_thickness(self):
        '''
        Returns the layer thickness (cm).
        '''
        sum = 0
        for layer in self:
            sum += layer.get_thickness()
        return sum

    def get_material(self, x):
        '''
        Returns the material at a given position (cm) within the layer.
        
        Parameters
        ----------
        x : float
            position (cm) between zero and the layer thickness
        '''
        if x < 0. or x > self.get_thickness():
            raise ValueError('x not within range [{:g}, {:g}]'
                             ''.format(0., self.get_thickness()))
        last_x = 0
        for layer in self:
            next_x = last_x + layer.get_thickness()
            if x <= next_x:
                return layer.get_material(x-last_x)
            last_x = next_x
        else:
            raise RuntimeError('unexpected execution path')

    def get_flatband(self, T=300.):
        '''
        returns x, Ev, Ec, Ei
        
        x will be [0, ..., thickness]
        Ev will be [VBO, ..., VBO]
        Ec will be [CBO, ..., CBO]
        Ei will be [VBO+Ei, ..., VBO+Ei]
        
        Arguments
        ---------
        T : float
            the temperature
        '''
        x, Ev, Ec, Ei = [], [], [], []
        last_x = 0
        for layer in self:
            l_x, l_Ev, l_Ec, l_Ei = layer.get_flatband(T=T)
            for layer_x in l_x:
                x.append(last_x + layer_x)
            Ev.extend(l_Ev)
            Ec.extend(l_Ec)
            Ei.extend(l_Ei)
            thickness = l_x[-1]
            last_x += thickness
        return x, Ev, Ec, Ei