from .units import cm
from .layer import Layer
from .contact import Contact
from .solver import poisson_eq


class LayerStructure(object):
    def __init__(self, items=[]):
        self._items = []
        for item in items:
            self.append(item)

    def __iter__(self):
        return self._items.__iter__()

    def insert(self, index, item):
        if not (isinstance(item, Layer) or
                isinstance(item, LayerStructure) or
                isinstance(item, Contact)):
            raise ValueError('unkown item type: {}'.format(type(item)))
        self._items.insert(index, item)

    def append(self, item):
        if not (isinstance(item, Layer) or
                isinstance(item, LayerStructure) or
                isinstance(item, Contact)):
            raise ValueError('unkown item type: {}'.format(type(item)))
        self._items.append(item)

    def get_thickness(self):
        '''
        Returns the thickness (cm)
        '''
        sum = 0
        for item in self._items:
            if not isinstance(item, Contact):
                sum += item.get_thickness()
        return sum

    def get_material_at_depth(self, x):
        last_x = 0
        for item in self._items:
            if not isinstance(item, Contact):
                new_x = last_x + item.get_thickness()
                if new_x >= x:
                    if not isinstance(item, Layer):
                        return item.get_material_at_depth(x - last_x)
                    else:
                        return item.material
                last_x = new_x

    def write_AMBER_recipe(self, fobj):
        '''Writes the AMBER recipe commands to the file object'''
        fobj.write('! ** Total thickness: {:.0f} Angs\n'
                   ''.format(self.get_thickness() * 1e8))
        for layer, repeat in reversed(zip(self.layers, self.repeats)):
            if repeat > 1:
                fobj.write('! ' + '*'*80 + '\n')
                fobj.write('! * Begin repeated structure, {:.0f} times\n'.format(repeat))
                fobj.write('! ' + '*'*80 + '\n')
                fobj.write('repeat {:.0f}\n'.format(repeat))
            layer.write_commands(fobj)
            if repeat > 1:
                fobj.write('er\n'.format(repeat))
                fobj.write('! ' + '*'*80 + '\n')
                fobj.write('! * End repeated structure,  {:.0f} times\n'.format(repeat))
                fobj.write('! ' + '*'*80 + '\n')

    def get_flatband(self, T=300):
        '''
        returns x, Ev, Ec
        
        x will be [0, ..., thickness]
        Ev will be [VBO, ..., VBO]
        Ec will be [VBO+Eg, ..., VBO+Eg]
        
        Arguments
        ---------
        T : float
            the temperature
        '''
        x, Ev, Ec, Ei = [], [], [], []
        _x = 0
        for item in self._items:
            if isinstance(item, Contact):
                continue
            else:
                l_x, l_Ev, l_Ec, l_Ei = item.get_flatband(T=T)
                x.extend([x_+_x for x_ in l_x])
                Ev.extend([E_ for E_ in l_Ev])
                Ec.extend([E_ for E_ in l_Ec])
                Ei.extend([E_ for E_ in l_Ei])
                thickness = item.get_thickness()
                _x += thickness
        return x, Ev, Ec, Ei

    def show_flatband(self, T=300):
        '''
        Show a plot of the band profile at flatband.
        
        Arguments
        ---------
        T : float
            the temperature
        '''
        import matplotlib.pyplot as plt
        import numpy
        fig, ax = plt.subplots()
        x, Ev, Ec, Ei = self.get_flatband(T=T)
        x = numpy.array(x)*1e7 # nm
        ax.plot(x, Ev, 'r-', label='$E_v$')
        ax.plot(x, Ec, 'b-', label='$E_c$')
        ax.plot(x, Ei, 'k:', label='$E_i$')
        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel('Depth (nm)')
        plt.show()

    def show_equilibrium(self, T=300, N=1000):
        '''
        Show a plot of the band profile at equilibrium.
        
        Arguments
        ---------
        T : float
            the temperature
        N : int
            the number of grid points
        '''
        import numpy
        x, Ev, Ec, Ei, p, n, Na, Nd = poisson_eq(self, T=T, N=N)
        x = x*1e7 # nm
        Ef = [0 for x_ in x]
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='col')
        ax1.plot(x, Ev, 'r-', label='$E_v$')
        ax1.plot(x, Ec, 'b-', label='$E_c$')
        ax1.plot(x, Ef, 'k--', label='$E_f$')
        ax1.plot(x, Ei, 'k:', label='$E_i$')
        ax1.set_ylabel('Energy (eV)')
        ax2.semilogy(x, Na, 'r-', label='$N_A$')
        ax2.semilogy(x, Nd, 'b-', label='$N_D$')
        ax2.semilogy(x, p, 'r--', label='$p$')
        ax2.semilogy(x, n, 'b--', label='$n$')
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
    