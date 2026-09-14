#!/usr/bin/env python3
"""
One-at-a-time control-input sweeps (diabetes parameterisation, inflamed
initial state [0.7, 0.15, 0.20], t = 200).

Quantifies how strongly each control input can move the composite
regenerative score across its full operating range:

    k_ROS     in [0, 0.8] : RS 0.12 -> 0.92  (range 0.80)  <- dominant
    alpha_AMPK in [0, 0.5]: RS 0.12 -> 0.57  (range 0.44)
    E_matrix  in [0, 0.3] : RS 0.12 -> 0.30  (range 0.18)

This replaces the earlier (unverifiable) Sobol S_total values with a
directly reproducible measure of control authority.
"""

import numpy as np
from ode_model import DEFAULT_PARAMS_DM, simulate, get_final_state, regscore


def sweep(name, values, init=(0.7, 0.15, 0.20)):
    scores = []
    for v in values:
        sol = simulate(DEFAULT_PARAMS_DM, {name: float(v)}, list(init))
        X, Y, Z = get_final_state(sol)
        scores.append(regscore(X, Y, Z))
    return np.array(scores)


if __name__ == '__main__':
    for name, vals in [('k_ROS', np.linspace(0, 0.8, 41)),
                       ('alpha_AMPK', np.linspace(0, 0.5, 41)),
                       ('E_matrix', np.linspace(0, 0.3, 41))]:
        rs = sweep(name, vals)
        print(f"{name:12s}: RS {rs.min():.3f} -> {rs.max():.3f} "
              f"(achievable range {rs.max() - rs.min():.3f})")
