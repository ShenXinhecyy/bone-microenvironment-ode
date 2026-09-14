#!/usr/bin/env python3
"""
Bifurcation analysis of the BOLD model (diabetes parameterisation).

Two complementary computations:
  1. Equilibrium continuation (Newton root-finding via fsolve) of all three
     steady-state branches across k_ROS in [0, 0.8]. Locates the saddle-node
     fold of the inflammatory branch at k_ROS* = 0.333.
  2. Finite-time grid scan (simulation to t = 200 from the inflamed initial
     state [0.7, 0.15, 0.20]), the numerically less precise but more
     intuitive estimator; locates the transition at k_ROS = 0.338.

Outputs: Figure4_bifurcation_data.csv (k_ROS, RS_inflammatory, RS_unstable,
RS_regenerative).
"""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from ode_model import DEFAULT_PARAMS_DM, ode_system, simulate, get_final_state, regscore


def steady(control, guess, params=DEFAULT_PARAMS_DM):
    """Newton solve for a steady state near `guess`; returns None if not physical."""
    x, info, ier, msg = fsolve(lambda s: ode_system(0, s, params, control),
                               guess, full_output=True)
    if ier == 1 and np.all(x > -1e-6) and np.all(x < 1 + 1e-6):
        return np.clip(x, 0, 1)
    return None


def equilibrium_bifurcation(path='Figure4_bifurcation_data.csv'):
    ks = np.round(np.arange(0, 0.801, 0.005), 3)
    rows = []
    for k in ks:
        c = {'k_ROS': float(k)}
        s_inf = steady(c, [0.79, 0.14, 0.06])   # inflammatory branch
        s_reg = steady(c, [0.05, 0.80, 0.95])   # regenerative branch
        s_mid = steady(c, [0.35, 0.35, 0.35])   # unstable (saddle) branch
        # keep only genuine branch members (Newton can jump basins)
        if s_inf is not None and not s_inf[0] > 0.4:
            s_inf = None
        if s_reg is not None and not (s_reg[2] > 0.8 and s_reg[0] < 0.2):
            s_reg = None
        rs_inf = regscore(*s_inf) if s_inf is not None else np.nan
        rs_reg = regscore(*s_reg) if s_reg is not None else np.nan
        rs_mid = regscore(*s_mid) if s_mid is not None else np.nan
        if s_mid is not None:
            if s_inf is not None and np.allclose(s_mid, s_inf, atol=1e-4):
                rs_mid = np.nan
            if s_reg is not None and np.allclose(s_mid, s_reg, atol=1e-4):
                rs_mid = np.nan
        rows.append((k, rs_inf, rs_mid, rs_reg))
    df = pd.DataFrame(rows, columns=['k_ROS', 'RS_inflammatory',
                                     'RS_unstable', 'RS_regenerative'])
    df.to_csv(path, index=False)
    fold = df.loc[df['RS_inflammatory'].notna(), 'k_ROS'].max() + 0.0025
    print(f"Equilibrium-continuation fold: k_ROS* = {fold:.4f}")
    return df


def grid_scan(init=(0.7, 0.15, 0.20), t_span=(0, 200)):
    """Finite-time grid scan (the 't = 200 grid' estimator)."""
    prev = None
    for k in np.arange(0, 0.801, 0.001):
        sol = simulate(DEFAULT_PARAMS_DM, {'k_ROS': float(k)}, list(init),
                       t_span=t_span)
        X, Y, Z = get_final_state(sol)
        rs = regscore(X, Y, Z)
        if rs > 0.5:
            print(f"Finite-time grid transition: k_ROS = {k:.3f}")
            return k
    return None


if __name__ == '__main__':
    equilibrium_bifurcation()
    grid_scan()
