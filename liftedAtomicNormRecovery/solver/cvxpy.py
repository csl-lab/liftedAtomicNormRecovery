## ----------------------------------------------------------------------------
## Imports

import numpy as np
import cvxpy as cp

from numpy.typing import NDArray
from scipy.sparse import spmatrix

CLARABEL_OPTS = { 
                  'max_iter'                : 1000000,
                  'tol_gap_abs'             :	1e-6,
                  'tol_gap_rel'	            : 1e-6,
                  'tol_feas'                :	1e-6,
                  'tol_infeas_abs'          : 1e-6,
                  'tol_infeas_rel'          : 1e-6,
                  'tol_ktratio'	            : 1e-4,
                  'reduced_tol_gap_abs'     : 5e-4,
                  'reduced_tol_gap_rel'	    : 5e-4,
                  'reduced_tol_feas'        : 1e-3,
                  'reduced_tol_infeas_abs'  : 5e-4,
                  'reduced_tol_infeas_rel'  : 5e-4,
                  'reduced_tol_ktratio'	    : 1e-3
                  }

OSQP_OPTS = { 
              'max_iter'                : 1000000
                  }

class AtomicNormRecovery(object):
  def __init__(self):
    self.threshold = 1E-6
    self.status = None
    self.verbose = False

  def solve(self, A : NDArray | spmatrix, T : NDArray | spmatrix, y : NDArray | spmatrix, delta : float = 0.0) -> NDArray:
    # clear
    self.status = 'Called solve.'
    # define variables
    cp_c = cp.Variable(A.shape[1], nonneg=True)
    # constraint matrix
    B = T @ A
    # define constraints
    if np.abs(delta) < self.threshold:
      cp_constraints = [ B @ cp_c == y ]
    else:
      cp_constraints = [ cp.norm(B @ cp_c - y, 2) <= delta ]
    # define objective
    cp_objective = cp.Minimize(cp.sum(cp_c))
    # define problem
    cp_problem = cp.Problem(cp_objective, cp_constraints)
    # solve
    cp_problem.solve(verbose=self.verbose, solver='CLARABEL', **CLARABEL_OPTS)
    # update
    self.status = cp_problem.solver_stats
    # retrieve optimal variable
    return cp_c.value
  
class LiftedAtomicNormRecoveryLiftCoordinate(object):
  def __init__(self):
    self.threshold = 1E-6
    self.status = None
    self.verbose = False

  def solve(self, A : NDArray | spmatrix, T : NDArray | spmatrix, y : NDArray | spmatrix, lam : NDArray | spmatrix, alpha : float = 0.0, delta : float = 0.0) -> NDArray:
    # clear
    self.status = 'Called solve.'
    if A.shape[0] == T.shape[1] + 1:
      lam = A[-1]
      A = A[:T.shape[1]]
    # define variables
    cp_c = cp.Variable(A.shape[1], nonneg=True)
    # cost function
    w = np.ones((A.shape[1],), dtype=float) - lam
    # constraint matrix
    B = T @ A
    # define constraints
    if np.abs(delta) < self.threshold:
      cp_constraints = [ B @ cp_c == y, lam @ cp_c <= alpha ]
    else:
      cp_constraints = [ cp.norm(B @ cp_c - y, 2) <= delta, lam @ cp_c <= alpha ]
    # define objective
    cp_objective = cp.Minimize(w @ cp_c)
    # define problem
    cp_problem = cp.Problem(cp_objective, cp_constraints)
    # solve
    try:
      cp_problem.solve(verbose=self.verbose, ignore_dpp=True, solver='CLARABEL', **CLARABEL_OPTS)
    except:
      # relax constraints
      cp_constraints = [ cp.norm(B @ cp_c - y, 2) <= 1E-5, lam @ cp_c <= alpha ]
      # define problem
      cp_problem = cp.Problem(cp_objective, cp_constraints)
      # solve
      cp_problem.solve(verbose=self.verbose, ignore_dpp=True, solver='CLARABEL', **CLARABEL_OPTS)
    # update
    self.status = cp_problem.solver_stats
    # retrieve optimal variable
    return cp_c.value