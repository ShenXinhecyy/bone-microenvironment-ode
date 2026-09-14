#!/usr/bin/env python3
"""
Core ODE model for bone regeneration microenvironment (BOLD framework).
Three variables: X (inflammation/ROS), Y (osteogenesis), Z (metabolic fitness;
M2-like macrophage + MSC metabolism). Unified saturation mechanism (Model C).

"""

import numpy as np
from scipy.integrate import solve_ivp


def hill(x, k, n):
    """Safe Hill function, avoids division by zero."""
    if x <= 0:
        return 0.0
    xn = x ** n
    kn = k ** n
    return xn / (kn + xn)


def regscore(X, Y, Z):
    """Composite regenerative score: weighted combination of the three variables."""
    return 0.4 * Z + 0.4 * Y + 0.2 * (1 - X)


def ode_system(t, state, params, control):
    """
    ODE system with unified saturation (Model C).

    state   : [X, Y, Z]
    params  : 20 model parameters (see DEFAULT_PARAMS_DM)
    control : external control inputs (k_ROS, alpha_AMPK, E_matrix)
    """
    X, Y, Z = state

    aX = params['alpha_X']
    bX = params['beta_X']
    gXY = params['gamma_XY']
    gZX = params['gamma_ZX']
    aYb = params['alpha_Y_base']
    gYZ = params['gamma_YZ']
    kXY = params['kappa_XY_inhib']
    dY_ = params['delta_Y']
    bY = params['beta_Y']
    aZ = params['alpha_Z']
    bZ = params['beta_Z']
    gZXd = params['gamma_ZX_damage']
    n, m, p, q = params['n'], params['m'], params['p'], params['q']
    KX, KXY, KZ, KX2 = params['K_X'], params['K_XY'], params['K_Z'], params['K_X2']

    kROS = control.get('k_ROS', 0.0)
    aAMPK = control.get('alpha_AMPK', 0.0)
    Emat = control.get('E_matrix', 0.0)

    H_X = hill(X, KX, n)
    H_Y = hill(Y, KXY, m)
    H_Z = hill(Z, KZ, q)
    H_X2 = hill(X, KX2, p)

    dX = aX * H_X * (1 - X) - bX * X - gXY * H_Y * X - gZX * H_Z * X - kROS * X
    dY = (aYb + gYZ * Z + bY * Y + Emat) * (1 - Y) - kXY * H_X2 * Y - dY_ * Y
    dZ = -gZXd * X * Z + (aZ * (1 - X) + bZ * Z + aAMPK) * (1 - Z)

    return [dX, dY, dZ]


# Default parameters for Diabetes Mellitus (DM) model — 20 parameters
# (12 rate/coupling constants + 4 Hill coefficients + 4 half-activation constants)
DEFAULT_PARAMS_DM = {
    'alpha_X': 0.80, 'beta_X': 0.15,
    'gamma_XY': 0.35, 'gamma_ZX': 0.30,
    'alpha_Y_base': 0.06, 'gamma_YZ': 0.25,
    'kappa_XY_inhib': 0.50,
    'delta_Y': 0.10, 'beta_Y': 0.10,
    'alpha_Z': 0.05, 'beta_Z': 0.08,
    'gamma_ZX_damage': 0.30,
    'n': 2.5, 'm': 2.0, 'p': 2.5, 'q': 2.0,
    'K_X': 0.35, 'K_XY': 0.45, 'K_Z': 0.42, 'K_X2': 0.35
}


def simulate(params, control, init_state, t_span=(0, 200), t_eval=None):
    """Run a single simulation (RK45, rtol=1e-8)."""
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 10)

    sol = solve_ivp(
        ode_system, t_span, init_state,
        args=(params, control),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        max_step=1.0
    )
    return sol


def get_final_state(sol):
    """Extract final (X, Y, Z) from solution."""
    return sol.y[0, -1], sol.y[1, -1], sol.y[2, -1]
