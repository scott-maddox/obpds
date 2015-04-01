from pint import UnitRegistry
units = UnitRegistry()

m = units.m
cm = units.cm
cm3 = units.cm**3
um = units.um
nm = units.nm

def to_units(v, u):
    if not isinstance(v, units.Quantity):
        raise TypeError('Missing units: {}'.format(v))
    return v.to(u).magnitude