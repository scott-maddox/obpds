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
from numpy.linalg import norm
from scipy.sparse.linalg import spsolve

class NewtonResult(object):
    def __init__(self, value, num_iter, converged):
        self.value = value
        self.num_iter = num_iter
        self.converged = converged

def scalar_newton(G, A, u0, rtol=0, atol=1e-6, tau=0.5, max_iter=100):
    '''
    Uses Newton's method to solve a non-linear equation,

        G(U) = 0,

    where U is the unkown scalar. The derivative,

        A(U) = dG(U)/dU,
    
    should be linearized about U. The Newton procedure used is described in
    the user manual.

    Arguments
    ---------
    G : callable
        Accepts the current guess, uk.
        Returns the residual scalar at uk.
    A : callable
        Accepts the current guess, uk.
        Returns the derivative linearized about uk.
    u0 : float
        Initial guess of the unkown scalar, U.
    rtol : float (default: 0)
        Relative tolerance for termination.
    atol : float (default: 1e-6)
        Absolute tolerance for termination.
    tau : float (default: 0.5)
        Decrease parameter, used in the algorithm to stabilize the Newton
        iteration. The value should be between 0 and 1. Values closer to 1
        increase stability, but slow down convergence.
    max_iter : int (default: 100)
        Maximum number of iterations.

    Returns
    -------
    result.value : float
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
    Gk_norm = numpy.abs(Gk)

    for k in xrange(1, max_iter+1):
        # Step 2. Calculate the delta.
        duk = -Gk/Ak
#         logger.debug('duk = %s', duk)

        aerr = numpy.abs(duk)
        logger.debug('aerr = %e', aerr)
        if aerr < atol:
            logger.info('converged after %i iterations', k)
            return NewtonResult(value=uk + duk, num_iter=k, converged=True)
        
        while True:
            # Step 3. Calculate the next guess.
            ukp = uk + sk*duk
            assert not numpy.isnan(ukp)

            Gkp = G(ukp)
            Gkp_norm = numpy.abs(Gkp)
            xi_kp = Gkp_norm / Gk_norm
            
            # Step 4. Adjust the damping parameter.
            if 1 - xi_kp < tau * sk:
                sk /= 2.
                logger.debug('decreasing sk to %e', sk)
                continue  # go to Step 3.
            else:
                sk = sk/(sk+(1.-sk)*xi_kp/2.)
                logger.debug('increasing sk to %e', sk)
                break  # go to Step 5.

        # Step 5. Quit if converged; iterate if not.
        rerr = numpy.abs(duk/ukp)
        logger.debug('rerr = %f', rerr)
        if rerr < rtol:
            logger.info('converged after %i iterations', k)
            return NewtonResult(value=ukp, num_iter=k, converged=True)
        else:
            uk = ukp
            Ak = A(uk)
            Gk = Gkp
            Gk_norm = Gkp_norm
            continue  # go to Step 2.

    # Failed to converge, but return what we have.
    logger.warn('failed to converge after %i iterations', k)
    return NewtonResult(value=ukp, num_iter=k, converged=False)


def newton(G, A, u0, rtol=0, atol=1e-6, tau=0.5, max_iter=100):
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
    rtol : float (default: 0)
        Relative tolerance for termination.
    atol : float (default: 1e-6)
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
#         logger.debug('duk = %s', duk)

        # Check for convergence using absolute tolerance
        aerr = numpy.abs(duk).max()
        logger.debug('aerr = %f', aerr)
        if aerr < atol:
            logger.info('converged after %i iterations', k)
            return NewtonResult(value=uk + duk, num_iter=k, converged=True)
        
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
                logger.debug('decreasing sk to %e', sk)
                continue  # go to Step 3.
            else:
                sk = sk/(sk+(1.-sk)*xi_kp/2.)
                logger.debug('increasing sk to %e', sk)
                break  # go to Step 5.

        # Step 5. Quit if converged; iterate if not.
        # Check for convergence using relative tolerance
        rerr = numpy.abs(duk/ukp).max()
        logger.debug('rerr = %e', rerr)
        if rerr < rtol:
            logger.info('converged after %i iterations', k)
            return NewtonResult(value=ukp, num_iter=k, converged=True)
        else:
            uk = ukp
            Ak = A(uk)
            Gk = Gkp
            Gk_norm = Gkp_norm
            continue  # go to Step 2.

    # Failed to converge, but return what we have.
    logger.warn('failed to converge after %i iterations', k)
    return NewtonResult(value=ukp, num_iter=k, converged=False)
