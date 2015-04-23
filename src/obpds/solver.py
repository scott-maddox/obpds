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

import numpy
from numpy import exp, log, sqrt
from numpy.linalg import norm
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix, spdiags

from .solution import EquilibriumSolution
from .contact import OhmicContact, SchottkyContact

class NewtonResult(object):
    def __init__(self, value, num_iter, converged):
        self.value = value
        self.num_iter = num_iter
        self.converged = converged

def newton(G, A, u0, atol=1e-4, tau=0.5, max_iter=100):
    '''
    Uses Newton's method to solve a system of non-linear equations,

        G(U) = 0,

    where the unkown vector, U, corresponds to the values of the finite
    element or volume solution, u_h, at the node points. The Jacobian matrix,

        A(U) = dG(U)/dU,
    
    is a sparse stiffness matrix corresponding to a linear elliptic
    boundary value problem (linearized about U). The Newton procedure used
    is described in the user manual.

    Arguments
    ---------
    G : callable
        Accepts the current guess, uk.
        Returns the residual vector of length N at uk.
    A : callable
        Accepts the current guess, uk.
        Returns the NxN Jacobian matrix linearized about uk.
    u0 : vector of length N
        Initial guess of the unkown vector, U.
    atol : float (default: 1e-4)
        Absolute tolerance for termination.
    tau : float (default: 0.5)
        Decrease parameter, used in the algorithm to stabilize the Newton
        iteration. The value should be between 0 and 1. Values closer to 1
        increase stability, but slow down convergence.
    max_iter : int (default: 100)
        Maximum number of iterations.

    Returns
    -------
    result.value : vector of length N
        The unkown vector, U.
    result.num_iter: int
        The number of iterations used.
    result.converged: 
        True if the solution converged, False if it did not.
    '''
    # Step 1. Calculate the initial values.
    sk = 1. # damping parameter
    uk = u0
    Ak = A(uk)
    Gk = G(uk)
    Gk_norm = norm(Gk)

    for k in xrange(1, max_iter+1):
        # Step 2. Calculate the delta.
        duk = spsolve(Ak, -Gk, use_umfpack=True)
#         print 'duk', duk

        while True:
            # Step 3. Calculate the next guess.
            ukp = uk + sk*duk
            assert not numpy.isnan(ukp).any()

            Gkp = G(ukp)
            Gkp_norm = norm(Gkp)
            xi_kp = Gkp_norm / Gk_norm
            
            # Step 4. Adjust the damping parameter.
            if 1 - xi_kp < tau * sk:
                sk /= 2.
#                 print 'decreasing sk to {}'.format(sk)
                continue  # go to Step 3.
            else:
                sk = sk/(sk+(1.-sk)*xi_kp/2.)
#                 print 'increasing sk to {}'.format(sk)
                break  # go to Step 5.

        # Step 5. Quit if converged; iterate if not.
        duk_max = numpy.abs(duk).max()
#         print duk_max, atol
        if duk_max < atol:
            print 'newton_poisson converged after {} iterations'.format(k)
            return NewtonResult(value=ukp, num_iter=k, converged=True)
        else:
            uk = ukp
            Ak = A(uk)
            Gk = Gkp
            Gk_norm = Gkp_norm
            continue  # go to Step 2.

    # Failed to converge, but return what we have.
    print 'WARNING: newton_poisson failed to converge after {} iterations'.format(k)
    return NewtonResult(value=ukp, num_iter=k, converged=False)


#TODO: shouldn't this depend on the gradient of eps_r?
def Fpsi(psi, f, dx, a, b):
    '''
    Used to solve the one-dimensional electrostatic Poisson equation,

        u'' = -f,

    under Dirichlet (fixed potential) boundary conditions:

        u[0] = a
        u[-1] = b

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

def jacobian__Fpsi__psi(df_du, dx):
    '''
    Used to solve the one-dimensional electrostatic Poisson equation,

        u'' = -f,

    under Dirichlet (fixed potential) boundary conditions:

        u[0] = a
        u[-1] = b

    Arguments
    ---------
    df_du : vector of length N
        derivative of the source term with respect to u [V cm**-3]
    dx : float
        grid spacing [cm]

    Returns
    -------
    jac__Fpsi__psi : square matrix of width/height N
        the Jacobian matrix
    '''
    N = df_du.size
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
        rows.append(i); cols.append( i ); data.append(-2 + df_du[i]*dx2)
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
    Vt = k*T  # eV
    flatband = device._get_flatband(T, N)
    x = flatband.x
    materials = device._get_materials(N)
    dx = x[1]  # cm
    
    # materials parameters
    eps_r = numpy.array([m.dielectric(T=T) for m in materials])
    eps = eps0 * eps_r  # C Vt**-1 cm**-1
    q_over_eps = q / eps  # Vt cm
    Na = numpy.array([m.Na(T=T) for m in materials])  # cm**-3
    Nd = numpy.array([m.Nd(T=T) for m in materials])  # cm**-3
    Nnet = numpy.array([m.Nnet(T=T) for m in materials])  # cm**-3
    Nc_Gamma = numpy.array([m.Nc_Gamma(T=T) for m in materials])  # cm**-3
    Nc_L = numpy.array([m.Nc_L(T=T) for m in materials])  # cm**-3
    Nc_X = numpy.array([m.Nc_X(T=T) for m in materials])  # cm**-3
    Nc = numpy.array([m.Nc(T=T) for m in materials])  # cm**-3
    Nv = numpy.array([m.Nv(T=T) for m in materials])  # cm**-3

    # band edges
    Ev0 = flatband.Ev  # cm**-3
    Ec_Gamma0 = flatband.Ec_Gamma  # cm**-3
    Ec_L0 = flatband.Ec_L  # cm**-3
    Ec_X0 = flatband.Ec_X  # cm**-3
    Ec0 = flatband.Ec  # cm**-3
    Ei0 = flatband.Ei
    
    # define functions for calculating p, n, dp, and dn
    if approx == 'boltzmann':
        def p(u):
            Ev = Ev0 - u
            return boltz_p(0., Ev, Nv, Vt)
        def n_Gamma(u):
            Ec_Gamma = Ec_Gamma0 - u
            return boltz_n(0., Ec_Gamma, Nc_Gamma, Vt)
        def n_L(u):
            Ec_L = Ec_L0 - u
            return boltz_n(0., Ec_L, Nc_L, Vt)
        def n_X(u):
            Ec_X = Ec_X0 - u
            return boltz_n(0., Ec_X, Nc_X, Vt)
        def n(u):
            return n_Gamma(u)+n_L(u)+n_X(u)
        def dp(u):
            Ev = Ev0 - u
            return dboltz_p(0., Ev, Nv, Vt)
        def dn(u):
            Ec_Gamma = Ec_Gamma0 - u
            Ec_L = Ec_L0 - u
            Ec_X = Ec_X0 - u
            return (dboltz_n(0., Ec_Gamma, Nc_Gamma, Vt) +
                    dboltz_n(0., Ec_L, Nc_L, Vt) + 
                    dboltz_n(0., Ec_X, Nc_X, Vt))
    elif approx == 'parabolic':
        def p(u):
            Ev = Ev0 - u
            return fermi_p(0., Ev, Nv, Vt)
        def n_Gamma(u):
            Ec_Gamma = Ec_Gamma0 - u
            return fermi_n(0., Ec_Gamma, Nc_Gamma, Vt)
        def n_L(u):
            Ec_L = Ec_L0 - u
            return fermi_n(0., Ec_L, Nc_L, Vt)
        def n_X(u):
            Ec_X = Ec_X0 - u
            return fermi_n(0., Ec_X, Nc_X, Vt)
        def n(u):
            return n_Gamma(u)+n_L(u)+n_X(u)
        def dp(u):
            Ev = Ev0 - u
            return dfermi_p(0., Ev, Nv, Vt)
        def dn(u):
            Ec_Gamma = Ec_Gamma0 - u
            Ec_L = Ec_L0 - u
            Ec_X = Ec_X0 - u
            return (dfermi_n(0., Ec_Gamma, Nc_Gamma, Vt) +
                    dfermi_n(0., Ec_L, Nc_L, Vt) + 
                    dfermi_n(0., Ec_X, Nc_X, Vt))
    elif approx == 'kane':
        alpha = numpy.array([m.nonparabolicity(T=T) for m in materials])
        def p(u):
            Ev = Ev0 - u
            return fermi_p(0., Ev, Nv, Vt)
        def n_Gamma(u):
            Ec_Gamma = Ec_Gamma0 - u
            return npfermi_n(0., Ec_Gamma, Nc_Gamma, alpha, Vt)
        def n_L(u):
            Ec_L = Ec_L0 - u
            return fermi_n(0., Ec_L, Nc_L, Vt)
        def n_X(u):
            Ec_X = Ec_X0 - u
            return fermi_n(0., Ec_X, Nc_X, Vt)
        def n(u):
            return n_Gamma(u)+n_L(u)+n_X(u)
        def dp(u):
            Ev = Ev0 - u
            return dfermi_p(0., Ev, Nv, Vt)
        def dn(u):
            Ec_Gamma = Ec_Gamma0 - u
            Ec_L = Ec_L0 - u
            Ec_X = Ec_X0 - u
            return (dnpfermi_n(0., Ec_Gamma, Nc_Gamma, alpha, Vt) +
                    dfermi_n(0., Ec_L, Nc_L, Vt) + 
                    dfermi_n(0., Ec_X, Nc_X, Vt))
    else:
        raise ValueError('Invalid value for approx: {}'.format(approx))

    # band edge potentials
#     Delta_Egc = 0.  # TODO: include bandgap reduction
#     Delta_Egv = 0.  # TODO: include bandgap reduction
    Vn = Vt*log(Nc/Nc[0]) - (Ec0 - Ec0[0])# + Delta_Egc
    Vp = Vt*log(Nv/Nv[0]) - (Ev0 - Ev0[0])# + Delta_Egv

    # effective intrinsic carrier concentration
    nieff = sqrt(Nc*Nv)*exp(-(Ec0-Ev0)/2./Vt)  # cm**-3
    
    # Use charge neutrality to guess psi
    psi0 = numpy.zeros(N)
    for i in xrange(N):
        if Nnet[i] <= 0.:
            p0 = -Nnet[i]/2. + sqrt((Nnet[i]/2.)**2+nieff[i]**2.)
            psi0[i] = -Vt*log(p0/Nv[i]) + Ev0[i]
        else:
            n0 = Nnet[i]/2. + sqrt((Nnet[i]/2.)**2+nieff[i]**2)
            psi0[i] = Vt*log(n0/Nc[i])+Ec0[i]
    # Refine with newton method, in case the single conduction band
    # approximation is inadequate
    def G0(u):
        return (p(u)-n(u)+Nnet)
    def A0(u):
        df_du = (dp(u) - dn(u))
        return spdiags(df_du, [0], N, N, format='csr')
    result = newton(G0, A0, psi0)
    psi0 = result.value  # eV
    
    # boundary conditions
    if isinstance(device._contacts[0], OhmicContact):
        a=psi0[0]
    elif isinstance(device._contacts[0], SchottkyContact):
        wf = device._contacts[0].work_function
        a = materials[0].electron_affinity(T=T) + Ec0[0] - wf
    else:
        raise RuntimeError('unexpected execution path')
    if isinstance(device._contacts[1], OhmicContact):
        b=psi0[-1]
    elif isinstance(device._contacts[1], SchottkyContact):
        wf = device._contacts[1].work_function
        b = materials[0].electron_affinity(T=T) + Ec0[0] - wf
    else:
        raise RuntimeError('unexpected execution path')

#     global last_u
#     last_u = psi0
    def G(u):
#         global last_u
#         plot(last_u, u, p(u), n(u), Nnet)
#         last_u = u
        f = (p(u)-n(u)+Nnet)*q_over_eps
        _Fpsi = Fpsi(u, f, dx, a=a, b=b)
        return _Fpsi
    def A(u):
        df_du = q_over_eps*(dp(u) - dn(u))
        return jacobian__Fpsi__psi(df_du, dx)
    
    result = newton(G, A, psi0)
    psi = result.value  # eV
    return EquilibriumSolution(T, N, x, Na, Nd,
                               Ev0 - psi, Ec_Gamma0 - psi, Ec_L0 - psi,
                               Ec_X0 - psi, Ec0-psi, Ei0-psi, psi,
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