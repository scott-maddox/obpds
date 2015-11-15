Open Band Parameters Device Simulator (OBPDS)
=============================================

A free, open-source technology computer aided design software for simulating
semiconductor structures and devices.

OBPDS currently provides 1D electrostatic simulation of
III-V compound semiconductor heterostructures, similar to Prof. William
Frensley's `Bandprof`_. Materials parameters are provided by the
`Open Band Parameters`_ sister project.

If you would like to try OBPDS before installing it,
check out the `interactive tutorial`_.

The `source code`_ is graciously hosted by GitHub.

.. _`BandProf`: https://courses.ece.ubc.ca/480/downloads.htm
.. _`Open Band Parameters`: http://github.com/scott-maddox/openbandparams
.. _`interactive tutorial`: http://mybinder.org/repo/scott-maddox/obpds-binder/tutorial.ipynb
.. _`source code`: http://github.com/scott-maddox/obpds

Installation
============

In order to use OBPDS, you must having a working `Python`_ distribution
installed. Python 3 support has not yet been tested, so Python 2.7 is
suggested.

.. _`Python`: https://www.python.org/download/

From PyPi
---------

This is the recommended method for installing OBPDS. `PyPi`_ is the python
package index, which contains many python packages that can be easily installed
with a single command. To install OBPDS from `PyPi`_, open up a command
prompt and run the following command::

    pip install fdint

.. _`PyPi`: http://pypi.python.org/pypi

From Source
-----------

In order to install OBPDS from source, you must fist install the
following prerequisite packages:

- Numpy_
- Scipy_
- Matplotlib_
- OpenBandParams_

.. _`Numpy`: http://docs.scipy.org/doc/numpy/user/install.html
.. _`Scipy`: http://www.scipy.org/install.html
.. _`Matplotlib`: http://matplotlib.org/users/installing.html
.. _`OpenBandParams`: http://scott-maddox.github.io/openbandparams/installation.html

Once these are installed, download the latest release ``.zip`` or ``.tar.gz``
source package from the `github releases`_ page, extract its contents, and run
``python setup.py install`` from within the extracted directory.

.. _`github releases`: http://github.com/scott-maddox/obpds/releases/latest

Documentation
=============

Interactive `documentation`_ is hosted through mybinder.org.
If you have difficulty accessing the interactive version,
a static version is available `here`_.

.. _`documentation`: http://mybinder.org/repo/scott-maddox/obpds-binder
.. _`here`: https://github.com/scott-maddox/obpds-binder/blob/master/index.ipynb
