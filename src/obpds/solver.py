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
from scipy.sparse import csr_matrix, spdiags

from .newton import scalar_newton, newton
from .solution import ZeroCurrentSolution, EquilibriumSolution
from .contact import OhmicContact, SchottkyContact

#TODO: shouldn't this depend on the gradient of eps_r?
def update_Fpsi(Fpsi, psi, f, dx, a, b):
    '''
    Used to solve the one-dimensional electrostatic Poisson equation,

        psi'' = -f,

    under Dirichlet (fixed potential) boundary conditions:

        psi[0] = a
        psi[-1] = b

    Arguments
    ---------
    Fpsi : vector of length N
        Vector to update.
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
    '''
    dx2 = dx*dx
    Fpsi[1:-1] = psi[:-2] - 2*psi[1:-1] + psi[2:] + f[1:-1]*dx2
    
    # Dirichlet boundary conditions
    Fpsi[0] = psi[0] - a
    Fpsi[-1] = psi[-1] - b

def get_Fpsi_jacobian(N):
    '''
    Returns a sparse matrix that can be used to solve the one-dimensional
    electrostatic Poisson equation,
    
        psi'' = -f,
    
    under Dirichlet (fixed potential) boundary conditions:

        psi[0] = a
        psi[-1] = b
    
    Before solving for the potential, the returned sparce matrix should be
    updated with `update_Fpsi_jacobian`.
    
    Parameters
    ----------
    N : int
        number of grid points

    Returns
    -------
    A : csr_matrix
        Jacobian matrix with width and height of N. This should be updated
        with `update_Fpsi_jacobian`, before solving for the potential.
    '''
    data = numpy.empty((N-2)*3+2)
    rows = numpy.empty((N-2)*3+2)
    cols = numpy.empty((N-2)*3+2)

    # Dirichlet boundary condition
    j = 0
    i = 0
    rows[i] = j; cols[i] = j; data[i] = 1; i += 1

    # The interior matrix elements
    for j in xrange(1, N-1):
        rows[i] = j; cols[i] = j-1; data[i] = 1; i += 1
        rows[i] = j; cols[i] =  j ; data[i] = 1; i += 1
        rows[i] = j; cols[i] = j+1; data[i] = 1; i += 1

    # Dirichlet boundary condition
    j = N-1
    rows[i] = j; cols[i] = j; data[i] = 1; i += 1
    
    A = csr_matrix((data, (rows, cols)), (N, N))
    return A

def update_Fpsi_jacobian(Fpsi_jacobian, df_dpsi, dx):
    '''
    Updates the sparse matrix from `get_Fpsi_jacobian`, for solving the
    one-dimensional electrostatic Poisson equation,

        psi'' = -f,

    under Dirichlet (fixed potential) boundary conditions:

        psi[0] = a
        psi[-1] = b

    Arguments
    ---------
    Fpsi_jacobian : csr_matrix
        Jacobian matrix to update (from `get_Fpsi_jacobian`).
    df_dpsi : vector of length N
        Derivative of the source term with respect to psi [V cm**-3].
    dx : float
        Grid spacing [cm].
    '''
    N = df_dpsi.size
    dx2 = dx*dx
    Fpsi_jacobian.data[2::3] = -2 + df_dpsi[1:-1]*dx2


from fdint.scfd import *

# electron charge
q = 1.602176565e-19 # C
# vacuum permittivity
eps0 = 8.854187817620e-14 # C V**-1 cm**-1
# Boltzmann constant
k = 8.6173324e-5 # eV K**-1

def get_fermi_functions(phi_p, phi_n, Ev, Nv,
                        Ec_Gamma, Ec_L, Ec_X,
                        Nc_Gamma, Nc_L, Nc_X,
                        nonparabolicity, Vt, approx='kane'):
    
    # define functions for calculating p, n, dp, and dn
    if approx == 'boltzmann':
        def p(psi, buf=None):
            return boltzmann_p(phi_p, psi, Ev, Nv, Vt, out=buf)
        def n_Gamma(psi, buf=None):
            return boltzmann_n(phi_n, psi, Ec_Gamma, Nc_Gamma, Vt, out=buf)
        def n_L(psi, buf=None):
            return boltzmann_n(phi_n, psi, Ec_L, Nc_L, Vt, out=buf)
        def n_X(psi, buf=None):
            return boltzmann_n(phi_n, psi, Ec_X, Nc_X, Vt, out=buf)
        def n(psi, buf=None):
            return boltzmann_n3(phi_n, psi, Ec_Gamma, Ec_L, Ec_X,
                             Nc_Gamma, Nc_L, Nc_X, Vt, out=buf)
        def dp(psi, buf=None):
            return boltzmann_dp(phi_p, psi, Ev, Nv, Vt, out=buf)
        def dn(psi, buf=None):
            return boltzmann_dn3(phi_n, psi, Ec_Gamma, Ec_L, Ec_X,
                             Nc_Gamma, Nc_L, Nc_X, Vt, out=buf)
    elif approx == 'parabolic':
        def p(psi, buf=None):
            return parabolic_p(phi_p, psi, Ev, Nv, Vt, out=buf)
        def n_Gamma(psi, buf=None):
            return parabolic_n(phi_n, psi, Ec_Gamma, Nc_Gamma, Vt, out=buf)
        def n_L(psi, buf=None):
            return parabolic_n(phi_n, psi, Ec_L, Nc_L, Vt, out=buf)
        def n_X(psi, buf=None):
            return parabolic_n(phi_n, psi, Ec_X, Nc_X, Vt, out=buf)
        def n(psi, buf=None):
            return parabolic_n3(phi_n, psi, Ec_Gamma, Ec_L, Ec_X,
                                Nc_Gamma, Nc_L, Nc_X, Vt, out=buf)
        def dp(psi, buf=None):
            return parabolic_dp(phi_p, psi, Ev, Nv, Vt, out=buf)
        def dn(psi, buf=None):
            return parabolic_dn3(phi_n, psi, Ec_Gamma, Ec_L, Ec_X,
                                 Nc_Gamma, Nc_L, Nc_X, Vt, out=buf)
    elif approx == 'kane':
        def p(psi, buf=None):
            return parabolic_p(phi_p, psi, Ev, Nv, Vt, out=buf)
        def n_Gamma(psi, buf=None):
            return nonparabolic_n(phi_n, psi, nonparabolicity,
                                  Ec_Gamma, Nc_Gamma, Vt, out=buf)
        def n_L(psi, buf=None):
            return parabolic_n(phi_n, psi, Ec_L, Nc_L, Vt, out=buf)
        def n_X(psi, buf=None):
            return parabolic_n(phi_n, psi, Ec_X, Nc_X, Vt, out=buf)
        def n(psi, buf=None):
            return nonparabolic_n3(phi_n, psi, nonparabolicity,
                                   Ec_Gamma, Ec_L, Ec_X,
                                   Nc_Gamma, Nc_L, Nc_X, Vt, out=buf)
        def dp(psi, buf=None):
            return parabolic_dp(phi_p, psi, Ev, Nv, Vt, out=buf)
        def dn(psi, buf=None):
            return nonparabolic_dn3(phi_n, psi, nonparabolicity,
                                    Ec_Gamma, Ec_L, Ec_X,
                                    Nc_Gamma, Nc_L, Nc_X, Vt, out=buf)
    else:
        raise ValueError('Invalid value for approx: {}'.format(approx))
    return p, n_Gamma, n_X, n_L, n, dp, dn

def scalar_charge_neutrality(psi0, p, n, dp, dn, Nnet):
    if Nnet == 0.:
        return psi0

    def G(psi):
        rho = p(psi)-n(psi)+Nnet
        # replace nan's with 0.0
        if numpy.isnan(rho):
            rho = 0.0
        logger.debug('rho = %s', rho)
        return rho
    def A(psi):
        drho_dpsi = dp(psi)-dn(psi)
        # replace nan's with 1.0
        if numpy.isnan(drho_dpsi) or drho_dpsi == 0.:
            drho_dpsi = 1.0
        logger.debug('drho_dpsi = %s', drho_dpsi)
        return drho_dpsi
    result = scalar_newton(G, A, psi0, rtol=1e-7)
    assert result.converged
    return result.value

def charge_neutrality(device, V, phi_p, phi_n, T=300., N=1000,
                      approx='kane'):
    # get materials and device parameters
    parameters = device._get_parameters(T, N)
    flatband = device._get_flatband(T, N)
    
    # calculate other parameters
    Vt = k*T  # eV
    
    # get functions for calculating p, n, dp, and dn
    p, _n_Gamma, _n_X, _n_L, n, dp, dn = get_fermi_functions(phi_p, phi_n,
                    flatband.Ev, parameters.Nv,
                    flatband.Ec_Gamma, flatband.Ec_L, flatband.Ec_X,
                    parameters.Nc_Gamma, parameters.Nc_L, parameters.Nc_X,
                    parameters.nonparabolicity, Vt, approx)

    ## Begin by setting the majority carrier density to the net doping density.
    psi0 = numpy.zeros(N)
    residual = numpy.zeros(N)
    for i in xrange(N):
        try:
            phi_pi = phi_p[i]
            phi_ni = phi_n[i]
        except:
            phi_pi = phi_p
            phi_ni = phi_n
        if parameters.Nnet[i] < 0.:
            p0 = -parameters.Nnet[i]
            phi_p0 = iparabolic_p(p0, 0., flatband.Ev[i], parameters.Nv[i], Vt)
            psi0[i] = phi_p0-phi_pi
        elif parameters.Nnet[i] > 0.:
            n0 = parameters.Nnet[i]
            phi_n0 = iparabolic_n(n0, 0., flatband.Ec[i], parameters.Nc[i], Vt)
            psi0[i] = phi_n0-phi_ni
        else:
            if phi_pi != numpy.inf and phi_ni != -numpy.inf:
                psi0[i] = flatband.Ei[i] - (phi_pi+phi_ni)/2
            elif phi_ni != -numpy.inf:
                psi0[i] = flatband.Ei[i] - phi_ni
                residual[i] = -parameters.ni[i]
            elif phi_pi != numpy.inf:
                psi0[i] = flatband.Ei[i] - phi_pi
                residual[i] = parameters.ni[i]
            else:
                psi0[i] = flatband.Ei[i]
    logger.debug('psi0 before newton = %s', psi0)
    ## Refine with newton method, in case satallite valleys or
    ## band non-parabolicity are significant.
    
    buf1 = numpy.empty(N)
    buf2 = numpy.empty(N)
#     global last_psi
#     last_psi = psi0
    def G0(psi):
#         global last_psi
#         plot(last_psi, psi, p(psi), n(psi), parameters.Nnet)
#         last_psi = psi
        rho = p(psi, buf1) - n(psi, buf2) + parameters.Nnet - residual
        # replace nan's with 0.0
        rho[numpy.isnan(rho)] = 0.0
        logger.debug('rho = %s', rho)
        return rho
    def A0(psi):
        drho_dpsi = dp(psi, buf1) - dn(psi, buf2)
        # replace nan's with 1.0
        drho_dpsi[numpy.isnan(drho_dpsi)] = 1.0
        logger.debug('df_dpsi = %s', drho_dpsi)
        return spdiags(drho_dpsi, [0], N, N, format='csr')
    result = newton(G0, A0, psi0)
    psi0 = result.value  # eV
    logger.debug('psi0 = %s', psi0)
    assert not numpy.isnan(psi0).any()
    return psi0

def _poisson_zero_current(device, psi0, phi_p, phi_n, V, T=300., N=1000,
                          approx='kane'):
    '''
    Uses Newton's method to solve the self-consistent electrostatic Poisson
    equation for the given device under equilibrium conditions.

    Arguments
    ---------
    device : TwoTerminalDevice
        Device
    psi0 : ndarray
        initial guess of the electrostatic potential
    phi_p : ndarray
        hole quasi-Fermi level at each grid point
    phi_n : ndarray
        electron quasi-Fermi level at each grid point
    T : float (default=300.)
        Device temperature
    N : int (default=1000)
        Number of grid points
    approx : str (default ='kane')
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
    p, n_Gamma, n_X, n_L, n, dp, dn = get_fermi_functions(phi_p, phi_n,
                    flatband.Ev, parameters.Nv,
                    flatband.Ec_Gamma, flatband.Ec_L, flatband.Ec_X,
                    parameters.Nc_Gamma, parameters.Nc_L, parameters.Nc_X,
                    parameters.nonparabolicity, Vt, approx)

    # band edge potentials
#     Delta_Egc = 0.  # TODO: include bandgap reduction
#     Delta_Egv = 0.  # TODO: include bandgap reduction
    
    # boundary conditions
    if isinstance(device._contacts[0], OhmicContact):
        p0, _, _, _, n0, dp0, dn0 = get_fermi_functions(phi_p[0], phi_n[0],
                flatband.Ev[0], parameters.Nv[0],
                flatband.Ec_Gamma[0], flatband.Ec_L[0], flatband.Ec_X[0],
                parameters.Nc_Gamma[0], parameters.Nc_L[0], parameters.Nc_X[0],
                parameters.nonparabolicity[0], Vt, approx)
        a = scalar_charge_neutrality(psi0[0], p0, n0, dp0, dn0,
                                     parameters.Nnet[0])
    elif isinstance(device._contacts[0], SchottkyContact):
        wf = device._contacts[0].work_function
        a = parameters.electron_affinity[0] + flatband.Ec[0] - wf
    else:
        raise RuntimeError('unexpected execution path')
    if isinstance(device._contacts[-1], OhmicContact):
        p0, _, _, _, n0, dp0, dn0 = get_fermi_functions(phi_p[-1], phi_n[-1],
                flatband.Ev[-1], parameters.Nv[-1],
                flatband.Ec_Gamma[-1], flatband.Ec_L[-1], flatband.Ec_X[-1],
                parameters.Nc_Gamma[-1], parameters.Nc_L[-1], parameters.Nc_X[-1],
                parameters.nonparabolicity[-1], Vt, approx)
        b = scalar_charge_neutrality(psi0[-1], p0, n0, dp0, dn0,
                                     parameters.Nnet[-1])
    elif isinstance(device._contacts[-1], SchottkyContact):
        wf = device._contacts[-1].work_function
        b = parameters.electron_affinity[-1] + flatband.Ec[-1] - wf #- V
    else:
        raise RuntimeError('unexpected execution path')

    buf1 = numpy.empty(N)
    buf2 = numpy.empty(N)
    Fpsi = numpy.empty(N)
    Fpsi_jacobian = get_Fpsi_jacobian(N)
#     global last_psi
#     last_psi = psi0
    def G(psi):
#         global last_psi
#         plot(last_psi, psi, p(psi), n(psi), parameters.Nnet)
#         last_psi = psi
        # f = q/eps * (p - n + Nnet)
        f = (p(psi, buf1) - n(psi, buf2) + parameters.Nnet) * q_over_eps
        # replace nan's with 0.0
        f[numpy.isnan(f)] = 0.0
        logger.debug('f = %s', f)
        update_Fpsi(Fpsi, psi, f, dx, a=a, b=b)
        return Fpsi
    def A(psi):
        # df_dpsi = q/eps * (dp - dn)
        df_dpsi = (dp(psi, buf1) - dn(psi, buf2)) * q_over_eps
        # replace nan's with 1.0
        df_dpsi[numpy.isnan(df_dpsi)] = 1.0
        logger.debug('df_dpsi = %s', df_dpsi)
        update_Fpsi_jacobian(Fpsi_jacobian, df_dpsi, dx)
        return Fpsi_jacobian
    
    result = newton(G, A, psi0)
    psi = result.value  # eV
    Ev = flatband.Ev - psi
    Ec_Gamma = flatband.Ec_Gamma - psi
    Ec_L = flatband.Ec_L - psi
    Ec_X = flatband.Ec_X - psi
    Ec = flatband.Ec-psi
    Ei = flatband.Ei-psi
    field = -numpy.gradient(psi)/dx
    dEv_dx = numpy.gradient(Ev)/dx
    dEc_Gamma_dx = numpy.gradient(Ec_Gamma)/dx
    dEc_L_dx = numpy.gradient(Ec_L)/dx
    dEc_X_dx = numpy.gradient(Ec_X)/dx
    dEc_dx = numpy.gradient(Ec)/dx
    logger.debug('psi = %s', psi)
    return ZeroCurrentSolution(V, T, N, x,
                               parameters.Na, parameters.Nd,
                               phi_p, phi_n,
                               Ev, Ec_Gamma, Ec_L, Ec_X, Ec, Ei,
                               dEv_dx, dEc_Gamma_dx, dEc_L_dx, dEc_X_dx, dEc_dx,
                               psi, field,
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

def poisson_eq(device, T=300., N=1000, approx='kane'):
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
    approx : str (default ='kane')
        If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
        bands approximation (fastest). If 'parabolic', use the parabolic
        bands approximation (fast). If 'kane', include Gamma-valley
        non-parabolicity under the k.p Kane approximation (slow).

    Returns
    -------
    s : EquilibriumSolution
        Equilibrium solution
    '''
    V = phi_p = phi_n = numpy.zeros(N)
    
    # Use charge neutrality to guess psi.
    psi0 = charge_neutrality(device, V, phi_p, phi_n, T, N, approx)
    
    zcs = _poisson_zero_current(device, psi0, phi_p, phi_n, 0., T, N, approx)
    return EquilibriumSolution(zcs)

def get_phis(device, V, T, N):
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
    
    return phi_p, phi_n

def poisson_zero_current(device, V, T=300., N=1000, approx='kane'):
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
    approx : str (default ='kane')
        If 'boltzmann', use the Boltzmann (non-degenerate) and parabolic
        bands approximation (fastest). If 'parabolic', use the parabolic
        bands approximation (fast). If 'kane', include Gamma-valley
        non-parabolicity under the k.p Kane approximation (slow).

    Returns
    -------
    s : ZeroCurrentSolution
        Zero current solution
    '''
    phi_p, phi_n = get_phis(device, V, T, N)
    
    # Use charge neutrality to guess psi.
    psi0 = charge_neutrality(device, V, phi_p, phi_n, T, N, approx)

    return _poisson_zero_current(device, psi0, phi_p, phi_n, V, T, N, approx)

def capacitance_zero_current(device, V, dV, T=300., N=1000, approx='kane'):
    '''
    Calculate the capacitance under the zero-current approximation.
    '''
    s1 = device.get_zero_current(V, T, N, approx)
    phi_p, phi_n = get_phis(device, V+dV, T, N)
    s2 = _poisson_zero_current(device, s1.psi, phi_p, phi_n, V+dV, T, N, approx)
    zeros = numpy.zeros(N)
    drho = ((s2.p - s2.n) - (s1.p - s1.n))*q
    # The capacitance should be calculated from the charge on only one plate.
    # We could take the absolute value, integrate, and then divide by 2, but
    # that doesn't handle one side of the junction extending to the edge of the
    # the simulated region. Instead, we take the larger of the positive or
    # negative dQ. This is necessary because one of the plates might reach
    # outside of the simulated region.
    drho_positive = numpy.where(drho > 0., drho, zeros)
    drho_negative = numpy.where(drho < 0., drho, zeros)
    dQ_positive = numpy.trapz(drho_positive, dx=s1.x[1])
    dQ_negative = numpy.trapz(drho_negative, dx=s1.x[1])
    dQ = max(dQ_positive, -dQ_negative)
    C = dQ / dV
    return C
