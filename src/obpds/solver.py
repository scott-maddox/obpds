import numpy
from numpy import exp, log, sqrt
from numpy.linalg import norm
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

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
def poisson_eq(ls, T=300, N=1000):
    '''
    Uses Newton's method to solve the self-consistent electrostatic Poisson
    equation for the given layer structure, ls, under equilibrium conditions
    at temperature, T.

    Arguments
    ---------
    ls : LayerStructure
        Layer structure
    T : float
        Temperature [K]
    N : int
        Number of uniformly spaced grid points

    Returns
    -------
    x : vector of length N
        Depth
    Ev : vector of length N
        Valance band energy [eV]
    Ec : vector of length N
        Conduction band energy [eV]
    Ei : vector of length N
        Intrinsic energy [eV]
    p : vector of length N
        Hole concentration [cm**-3]
    n : vector of length N
        Electron concentration [cm**-3]
    Na : vector of length N
        Ionized acceptor concentration [cm**-3]
    Nd : vector of length N
        Ionized donor concentration [cm**-3]
    '''
    Vt = k*T  # eV
    x = numpy.linspace(0., ls.get_thickness(), N)
    dx = x[1]  # cm
    mats = [ls.get_material_at_depth(x_i) for x_i in x]
    eps_r = numpy.array([mat.dielectric(T=T) for mat in mats])
    eps = eps0 * eps_r  # C Vt**-1 cm**-1
    q_over_eps = q / eps  # Vt cm
    Nnet = numpy.array([mat.Nnet(T=T) for mat in mats])  # cm**-3
    meff_e_DOS = numpy.array([mat.meff_e_DOS(T=T) for mat in mats])  # cm**-3
#     print 'meff_e_DOS', meff_e_DOS
    Nc = numpy.array([mat.Nc(T=T) for mat in mats])  # cm**-3
#     print 'Nc', Nc
    Ncref = Nc[0]  # cm**-3
    Nv = numpy.array([mat.Nv(T=T) for mat in mats])  # cm**-3
#     print 'Nv', Nv
    Nvref = Nv[0]  # cm**-3
    VBO = numpy.array([mat.VBO(T=T) for mat in mats])  # cm**-3
    VBOref = VBO[0]
    Eg = numpy.array([mat.Eg(T=T) for mat in mats])  # cm**-3
    Egref = Eg[0]
    Delta_Egc = 0  # TODO: include bandgap reduction
    Delta_Egv = 0  # TODO: include bandgap reduction
    Vn = Vt*log(Nc/Ncref) - (VBO-VBOref+Eg-Egref) + Delta_Egc
#     print 'Vn', Vn
    Vp = Vt*log(Nv/Nvref) + (VBO-VBOref) + Delta_Egv
#     print 'Vp', Vp
    nieff = sqrt(Nc*Nv)*exp(-(Eg-Vn-Vp)/2/Vt)  # cm**-3
    nieffref = nieff[0]  # cm**-3
    
    # Use charge neutrality to guess psi
    psi0 = numpy.zeros(N)
    for i in xrange(N):
        if Nnet[i] <= 0:
            p0 = -Nnet[i]/2 + sqrt((Nnet[i]/2)**2+nieff[i]**2)
            psi0[i] = -inv_fermi_p(psi=0, Vp=Vp[i], p=p0, nieffref=nieffref, Vt=Vt)
        else:
            n0 = Nnet[i]/2 + sqrt((Nnet[i]/2)**2+nieff[i]**2)
            psi0[i] = -inv_fermi_n(psi=0, Vn=Vn[i], n=n0, nieffref=nieffref, Vt=Vt)
#     print 'psi0', psi0
    a=psi0[0]
    b=psi0[-1]

    zero = numpy.zeros(N)
    global last_u 
    last_u = psi0
    def G(u):
        global last_u
        p = fermi_p(psi=u, Vp=Vp, phi_p=zero, nieffref=nieffref, Vt=Vt)
        n = fermi_n(psi=u, Vn=Vn, phi_n=zero, nieffref=nieffref, Vt=Vt)
#         plot(last_u, u, p, n, Nnet)
        f = (p-n+Nnet)*q_over_eps
        _Fpsi = Fpsi(u, f, dx, a=a, b=b)
        last_u = u
        return _Fpsi
    def A(u):
        p = fermi_p(psi=u, Vp=Vp, phi_p=zero, nieffref=nieffref, Vt=Vt)
        n = fermi_n(psi=u, Vn=Vn, phi_n=zero, nieffref=nieffref, Vt=Vt)
        df_du = -1.*(p+n)*q_over_eps/Vt
        return jacobian__Fpsi__psi(df_du, dx)
    
    result = newton(G, A, psi0)
    psi = result.value  # eV
    p = fermi_p(psi=psi, Vp=Vp, phi_p=zero, nieffref=nieffref, Vt=Vt)
    n = fermi_n(psi=psi, Vn=Vn, phi_n=zero, nieffref=nieffref, Vt=Vt)
    Ec = -Vt*log(n/Nc)
    Ev = Vt*log(p/Nv)
    Ei = numpy.array([Ev[i]+mats[i].Ei(T=T) for i in xrange(N)])
    Na = numpy.array([mats[i].Na(T=T) for i in xrange(N)])
    Nd = numpy.array([mats[i].Nd(T=T) for i in xrange(N)])
    return x, Ev, Ec, Ei, p, n, Na, Nd

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