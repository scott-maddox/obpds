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

from setuptools import setup, find_packages

# read in __version__
exec(open('src/obpds/version.py').read())

setup(
      name='obpds',
      version=__version__,  # read from version.py
      description='free, open-source technology computer aided design software'
                  ' for simulating semiconductor structures and devices',
      long_description=open('README.rst').read(),
      url='http://scott-maddox.github.io/obpds',
      author='Scott J. Maddox',
      author_email='smaddox@utexas.edu',
      license='AGPLv3',
      packages=['obpds',
                'obpds.tests',
                'obpds.examples'],
      package_dir={'obpds': 'src/obpds'},
      test_suite='nose.collector',
      install_requires=['fdint >= 2.0',
                        'openbandparams >= 0.9',
                        'scipy',
                        'matplotlib',
                        'pint',
                        'numpy', # install numpy first (for some reason this works)
                        ],
      )
