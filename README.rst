Open Band Parameters Device Simulator (OBPDS)
=============================================

A free, open-source technology computer aided design software for simulation
semiconductor structures and devices.

The plan is to begin by providing 1D, zero-current electrostatic simulation of
III-V compound semiconductor heterostructures, similar to Prof. William
Frensley's `Bandprof`_. Next, drift-diffusion current simulation will be added
under the assumption of small band discontinuities. Then large band
discontinuities will be accounted for with thermionic emission boundary
conditions.

Materials parameters are provided by the `Open Band Parameters`_ sister
project.

The `source code`_ and `documentation`_ (coming soon) are graciously hosted
by GitHub.

.. _`BandProf`: https://courses.ece.ubc.ca/480/downloads.htm
.. _`Open Band Parameters`: http://github.com/scott-maddox/openbandparams
.. _`source code`: http://github.com/scott-maddox/obpds
.. _`documentation`: http://scott-maddox.github.io/obpds


Installation
============

In order to use OBPDS, you must having a working `Python`_ distribution
installed. Python 3 support has not yet been tested, so Python 2.7 is
suggested.

.. _`Python`: https://www.python.org/download/

From PyPi
---------

This is the easiest method. Install from `PyPi`_ by running `pip install obpds`
from the command line.

.. _`PyPi`: http://pypi.python.org/pypi

From Github
-----------

First, you will need to install the following prerequisite packages:
<<<<<<< HEAD

=======
>>>>>>> 69b5bea992e066f48d2b205676aabd9b11a85010
- Numpy_
- Scipy_
- Matplotlib_
- OpenBandParams_

.. _`Numpy`: http://docs.scipy.org/doc/numpy/user/install.html
.. _`Scipy`: http://docs.scipy.org/doc/scipy/user/install.html
.. _`Matplotlib`: http://matplotlib.org/users/installing.html
.. _`OpenBandParams`: http://scott-maddox.github.io/openbandparams/installation.html

Once these are installed, download the latest release `.zip` or `.tar.gz`
source package from the `github page`_, extract its contents, and run
`python setup.py install` from within the extracted directory
(OBPDS is a pure-python library, for now, so no compiling occurs
during this installation).

.. _`github page`: http://github.com/scott-maddox/openbandparams/releases/latest

Documentation
=============

Once you have `openbandparams` installed, check out the `tutorial`_ to
get acquainted.

.. _`tutorial`: http://scott-maddox.github.io/obpds/tutorial