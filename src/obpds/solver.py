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

import logging
logger = logging.getLogger(__name__)

import numpy
from numpy import exp, log, sqrt
from scipy.sparse import csr_matrix, spdiags

from .newton import newton
from .solution import ZeroCurrentSolution, EquilibriumSolution
from .contact import OhmicContact, SchottkyContact

#TODO: shouldn't this depend on the gradient of eps_r?
def Fpsi(psi, f, dx, a, b):
    '''
    Used to solve the one-dimensional electrostatic Poisson equation,

        psi'' = -f,

    under Dirichlet (fixed potential) boundary conditions:

        psi[0] = a
        psi[-1] = b

    Arguments
    ---------
    psi : vector of length N
        last estimate of the electrostatic potential [V]
    f : vector of length N
        source term [V cm**-2]
    dx : scalar
        grid spacing [cm]
    a : float
        left Dirichlet boundary value [V]
    b : float
        right Dirichlet boundary value [V]

    Returns
    -------
    Fpsi : vector of length N
        residual electrostatic potential [V]
    '''
    dx2 = dx**2
    B = numpy.empty_like(psi)
    B[1:-1] = psi[:-2] - 2*psi[1:-1] + psi[2:] + f[1:-1]*dx2
    
    # Dirichlet boundary conditions
    B[0] = psi[0] - a
    B[-1] = psi[-1] - b
    return B

def jacobian__Fpsi__psi(df_dpsi, dx):
    '''
    Used to solve the one-dimensional electrostatic Poisson equation,

        psi'' = -f,

    under Dirichlet (fixed potential) boundary conditions:

        psi[0] = a
        psi[-1] = b

    Arguments
    ---------
    df_dpsi : vector of length N
        derivative of the source term with respect to psi [V cm**-3]
    dx : float
        grid spacing [cm]

    Returns
    -------
    jac__Fpsi__psi : square matrix of width/height N
        the Jacobian matrix
    '''
    N = df_dpsi.size
    dx2 = dx**2
    data = []
    rows = []
    cols = []

    # Dirichlet boundary condition
    i = 0
    rows.append(i); cols.append(i); data.append(1)

    # The interior matrix elements
    for i in xrange(1, N-1):
        rows.append(i); cols.append(i-1); data.append(1)
        rows.append(i); cols.append( i ); data.append(-2 + df_dpsi[i]*dx2)
        rows.append(i); cols.append(i+1); data.append(1)

    # Dirichlet boundary condition
    i = N-1
    rows.append(i); cols.append(i); data.append(1)
    
    A = csr_matrix((data, (rows, cols)), (N, N))
    return A


from .fermi import *

# electron charge
q = 1.602176565e-19 # C
# vacuum permittivity
eps0 = 8.854187817620e-14 # C V**-1 cm**-1
# Boltzmann constant
k = 8.6173324e-5 # eV K**-1

def get_fermi_functions(device, phi_p, phi_n, T=300., N=1000,
                        approx='parabolic'):
    # get materials and device parameters
    parameters = device._get_parameters(T, N)
    flatband = device._get_flatband(T, N)
    
    # calculate other parameters
    Vt = k*T  # eV
    
    # define functions for calculating p, n, dp, and dn
    if approx == 'boltzmann':
        def p(psi):
            Ev = flatband.Ev - psi
            return boltz_p(phi_p, Ev, parameters.Nv, Vt)
        def n_Gamma(psi):
            Ec_Gamma = flatband.Ec_Gamma - psi
            return boltz_n(phi_n, Ec_Gamma, parameters.Nc_Gamma, Vt)
        def n_L(psi):
            Ec_L = flatband.Ec_L - psi
            return boltz_n(phi_n, Ec_L, parameters.Nc_L, Vt)
        def n_X(psi):
            Ec_X = flatband.Ec_X - psi
            return boltz_n(phi_n, Ec_X, parameters.Nc_X, Vt)
        def n(psi):
            return n_Gamma(psi)+n_L(psi)+n_X(psi)
        def dp(psi):
            Ev = flatband.Ev - psi
            return dboltz_p(phi_p, Ev, parameters.Nv, Vt)
        def dn(psi):
            Ec_Gamma = flatband.Ec_Gamma - psi
            Ec_L = flatband.Ec_L - psi
            Ec_X = flatband.Ec_X - psi
            return (dboltz_n(phi_n, Ec_Gamma, parameters.Nc_Gamma, Vt) +
                    dboltz_n(phi_n, Ec_L, parameters.Nc_L, Vt) + 
                    dboltz_n(phi_n, Ec_X, parameters.Nc_X, Vt))
    elif approx == 'parabolic':
        def p(psi):
            Ev = flatband.Ev - psi
            return fermi_p(phi_p, Ev, parameters.Nv, Vt)
        def n_Gamma(psi):
            Ec_Gamma = flatband.Ec_Gamma - psi
            return fermi_n(phi_n, Ec_Gamma, parameters.Nc_Gamma, Vt)
        def n_L(psi):
            Ec_L = flatband.Ec_L - psi
            return fermi_n(phi_n, Ec_L, parameters.Nc_L, Vt)
        def n_X(psi):
            Ec_X = flatband.Ec_X - psi
            return fermi_n(phi_n, Ec_X, parameters.Nc_X, Vt)
        def n(psi):
            return n_Gamma(psi)+n_L(psi)+n_X(psi)
        def dp(psi):
            Ev = flatband.Ev - psi
            return dfermi_p(phi_p, Ev, parameters.Nv, Vt)
        def dn(psi):
            Ec_Gamma = flatband.Ec_Gamma - psi
            Ec_L = flatband.Ec_L - psi
            Ec_X = flatband.Ec_X - psi
            return (dfermi_n(phi_n, Ec_Gamma, parameters.Nc_Gamma, Vt) +
                    dfermi_n(phi_n, Ec_L, parameters.Nc_L, Vt) + 
                    dfermi_n(phi_n, Ec_X, parameters.Nc_X, Vt))
    elif approx == 'kane':
        def p(psi):
            Ev = flatband.Ev - psi
            return fermi_p(phi_p, Ev, parameters.Nv, Vt)
        def n_Gamma(psi):
            Ec_Gamma = flatband.Ec_Gamma - psi
            return npfermi_n(phi_n, Ec_Gamma, parameters.Nc_Gamma,
                             parameters.nonparabolicity, Vt)
        def n_L(psi):
            Ec_L = flatband.Ec_L - psi
            return fermi_n(phi_n, Ec_L, parameters.Nc_L, Vt)
        def n_X(psi):
            Ec_X = flatband.Ec_X - psi
            return fermi_n(phi_n, Ec_X, parameters.Nc_X, Vt)
        def n(psi):
            return n_Gamma(psi)+n_L(psi)+n_X(psi)
        def dp(psi):
            Ev = flatband.Ev - psi
            return dfermi_p(phi_p, Ev, parameters.Nv, Vt)
        def dn(psi):
            Ec_Gamma = flatband.Ec_Gamma - psi
            Ec_L = flatband.Ec_L - psi
            Ec_X = flatband.Ec_X - psi
            return (dnpfermi_n(phi_n, Ec_Gamma, parameters.Nc_Gamma,
                               parameters.nonparabolicity, Vt) +
                    dfermi_n(phi_n, Ec_L, parameters.Nc_L, Vt) + 
                    dfermi_n(phi_n, Ec_X, parameters.Nc_X, Vt))
    else:
        raise ValueError('Invalid value for approx: {}'.format(approx))
    return p, n_Gamma, n_X, n_L, n, dp, dn

def charge_neutrality(device, V, phi_p, phi_n, T=300., N=1000,
                      approx='parabolic'):
    # get materials and device parameters
    parameters = device._get_parameters(T, N)
    flatband = device._get_flatband(T, N)
    
    # calculate other parameters
    Vt = k*T  # eV

    ## Begin by setting the majority carrier density to the net doping density.
    psi0 = numpy.zeros(N)
    for i in xrange(N):
        if parameters.Nnet[i] < 0.:
            p0 = -parameters.Nnet[i]
            phi_p0 = ifermi_p(p0, flatband.Ev[i], parameters.Nv[i], Vt)
            try:
                psi0[i] = phi_p0-phi_p[i]
            except:
                psi0[i] = phi_p0-phi_p
        elif parameters.Nnet[i] > 0.:
            n0 = parameters.Nnet[i]
            phi_n0 = ifermi_n(n0, flatband.Ec[i], parameters.Nc[i], Vt)
            try:
                psi0[i] = phi_n0-phi_n[i]
            except:
                psi0[i] = phi_n0-phi_n
        else:
            psi0[i] = flatband.Ei[i]
    logger.debug('psi0 before newton = %s', str(psi0))
    ## Refine with newton method, in case satallite valleys or
    ## band non-parabolicity are significant.
    ### get functions for calculating p, n, dp, and dn
    p, n_Gamma, n_X, n_L, n, dp, dn = get_fermi_functions(device,
                                                          phi_p, phi_n,
                                                          T, N,
                                                          approx)
#     global last_psi
#     last_psi = psi0
    def G0(psi):
#         global last_psi
#         plot(last_psi, psi, p(psi), n(psi), parameters.Nnet)
#         last_psi = psi
        rho = p(psi)-n(psi)+parameters.Nnet
        return rho
    def A0(psi):
        df_dpsi = dp(psi) - dn(psi)
        return spdiags(df_dpsi, [0], N, N, format='csr')
    result = newton(G0, A0, psi0)
    psi0 = result.value  # eV
    logger.debug('psi0 = %s', str(psi0))
    assert not numpy.isnan(psi0).any()
    return psi0

def _poisson_zero_current(device, V, phi_p, phi_n, T=300., N=1000,
                          approx='parabolic'):
    '''
    Uses Newton's method to solve the self-consistent electrostatic Poisson
    equation for the given device under equilibrium conditions.

    Arguments
    ---------
    device : TwoTerminalDevice
        Device
    phi_p : float or ndarray
        hole quasi-Fermi level at each grid point
    phi_n : float or ndarray
        electron quasi-Fermi level at each grid point
    T : float (default=300.)
        Device temperature
    N : int (default=1000)
        Number of grid points
    approx : str (default='parabolic')
        If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
        bands approximation (fastest). If 'parabolic', use the parabolic
        bands approximation (fast). If 'kane', include Gamma-valley
        non-parabolicity under the k.p Kane approximation (slow).

    Returns
    -------
    s : ZeroCurrentSolution
        Zero current solution
    '''
    # get materials and device parameters
    parameters = device._get_parameters(T, N)
    flatband = device._get_flatband(T, N)
    
    # calculate other parameters
    x = flatband.x
    dx = x[1]  # cm
    Vt = k*T  # eV
    eps = eps0 * parameters.dielectric  # C Vt**-1 cm**-1
    q_over_eps = q / eps  # Vt cm
    
    # get functions for calculating p, n, dp, and dn
    p, n_Gamma, n_X, n_L, n, dp, dn = get_fermi_functions(device,
                                                          phi_p, phi_n,
                                                          T, N,
                                                          approx)

    # band edge potentials
#     Delta_Egc = 0.  # TODO: include bandgap reduction
#     Delta_Egv = 0.  # TODO: include bandgap reduction
    
    # Use charge neutrality to guess psi.
    psi0 = charge_neutrality(device, V, phi_p, phi_n, T, N, approx)
    
    # boundary conditions
    if isinstance(device._contacts[0], OhmicContact):
        a = psi0[0]
    elif isinstance(device._contacts[0], SchottkyContact):
        wf = device._contacts[0].work_function
        a = parameters.electron_affinity[0] + flatband.Ec[0] - wf
    else:
        raise RuntimeError('unexpected execution path')
    if isinstance(device._contacts[-1], OhmicContact):
        b = psi0[-1] #- V
    elif isinstance(device._contacts[-1], SchottkyContact):
        wf = device._contacts[-1].work_function
        b = parameters.electron_affinity[-1] + flatband.Ec[-1] - wf #- V
    else:
        raise RuntimeError('unexpected execution path')

#     global last_psi
#     last_psi = psi0
    def G(psi):
#         global last_psi
#         plot(last_psi, psi, p(psi), n(psi), parameters.Nnet)
#         last_psi = psi
        f = (p(psi)-n(psi)+parameters.Nnet)*q_over_eps
        _Fpsi = Fpsi(psi, f, dx, a=a, b=b)
        return _Fpsi
    def A(psi):
        df_dpsi = (dp(psi) - dn(psi))*q_over_eps
        return jacobian__Fpsi__psi(df_dpsi, dx)
    
    result = newton(G, A, psi0)
    psi = result.value  # eV
    logger.debug('psi = %s', str(psi))
    return ZeroCurrentSolution(T, N, x,
                               parameters.Na, parameters.Nd,
                               phi_p, phi_n,
                               flatband.Ev - psi,
                               flatband.Ec_Gamma - psi,
                               flatband.Ec_L - psi,
                               flatband.Ec_X - psi,
                               flatband.Ec-psi,
                               flatband.Ei-psi,
                               psi,
                               n_Gamma(psi), n_L(psi), n_X(psi), n(psi), p(psi))

def plot(psi, psi_kp, p, n, Nnet):
    import matplotlib.pyplot as plt
    x = numpy.arange(psi.size)
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.plot(x, -psi, 'k--', label='$E_i$')
    ax1.plot(x, -psi_kp, 'k:', label='$E_i^k+1$')
    ax1.legend(loc='right')
    ax2.semilogy(x, numpy.abs(Nnet), 'k-', label='$N_net$')
    ax2.semilogy(x, p, 'r--', label='$p$')
    ax2.semilogy(x, n, 'b--', label='$n$')
    ax2.legend(loc='right')
    plt.show()

def poisson_eq(device, T=300., N=1000, approx='parabolic'):
    '''
    Uses Newton's method to solve the self-consistent electrostatic Poisson
    equation for the given device under equilibrium conditions.

    Arguments
    ---------
    device : TwoTerminalDevice
        Device
    T : float (default=300.)
        Device temperature
    N : int (default=1000)
        Number of grid points
    approx : str (default='parabolic')
        If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
        bands approximation (fastest). If 'parabolic', use the parabolic
        bands approximation (fast). If 'kane', include Gamma-valley
        non-parabolicity under the k.p Kane approximation (slow).

    Returns
    -------
    s : EquilibriumSolution
        Equilibrium solution
    '''
    V = phi_p = phi_n = 0.
    zcs = _poisson_zero_current(device, V, phi_p, phi_n, T, N, approx)
    return EquilibriumSolution(zcs)

def poisson_zero_current(device, V, T=300., N=1000, approx='parabolic'):
    '''
    Uses Newton's method to solve the self-consistent electrostatic Poisson
    equation for the given device at a given bias voltage under the
    zero-current approximation.

    Arguments
    ---------
    device : TwoTerminalDevice
        Device
    V : float
        Bias voltage, i.e. left/top contact bias - right/bottom contact bias
    T : float (default=300.)
        Device temperature
    N : int (default=1000)
        Number of grid points
    approx : str (default='parabolic')
        If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
        bands approximation (fastest). If 'parabolic', use the parabolic
        bands approximation (fast). If 'kane', include Gamma-valley
        non-parabolicity under the k.p Kane approximation (slow).

    Returns
    -------
    s : EquilibriumSolution
        Equilibrium solution
    '''
    Vt = k*T  # eV
    phi_p = numpy.empty(N)
    phi_n = numpy.empty(N)

    Fp, Fn = device._get_Fp_Fn(N)
    for i in xrange(N):
        if Fp[i] is None:
            phi_p[i] = numpy.inf
        elif Fp[i] == 'left':
            phi_p[i] = 0.
        elif Fp[i] == 'right':
            phi_p[i] = V
        else:
            raise RuntimeError('Invalid value for Fp: {}'.format(Fp[i]))

        if Fn[i] is None:
            phi_n[i] = -numpy.inf
        elif Fn[i] == 'left':
            phi_n[i] = 0.
        elif Fn[i] == 'right':
            phi_n[i] = V
        else:
            raise RuntimeError('Invalid value for Fp: {}'.format(Fp[i]))

    if (phi_p == numpy.inf).all() and (phi_n == -numpy.inf).all():
        raise ValueError('Fp or Fn must be specified for at least one layer')

    return _poisson_zero_current(device, V, phi_p, phi_n, T, N, approx)