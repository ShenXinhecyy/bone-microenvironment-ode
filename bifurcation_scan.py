#!/usr/bin/env python3
"""
Bifurcation analysis: scan k_ROS (ROS clearance rate) as the control parameter
and observe the transition from high-inflammation to healing steady state.
"""

import numpy as np
from ode_model import DEFAULT_PARAMS_DM, simulate, get_final_state, regscore


def main():
    params = DEFAULT_PARAMS_DM
    init_state = [0.7, 0.15, 0.20]  # Start from inflamed state

    print("=" * 70)
    print("Bifurcation Scan — k_ROS as Control Parameter")
    print("=" * 70)

    k_range = np.linspace(0, 0.8, 17)
    print(f"{'k_ROS':>8s}  {'X':>8s}  {'Y':>8s}  {'Z':>8s}  {'RS':>8s}")
    print("-" * 45)

    prev_rs = None
    for k in k_range:
        control = {'k_ROS': k, 'alpha_AMPK': 0.0, 'E_matrix': 0.0}
        sol = simulate(params, control, init_state)
        X, Y, Z = get_final_state(sol)
        rs = regscore(X, Y, Z)

        marker = " <<<" if prev_rs is not None and abs(rs - prev_rs) > 0.3 else ""
        print(f"  {k:>6.2f}   {X:>7.4f}  {Y:>7.4f}  {Z:>7.4f}  {rs:>7.4f}{marker}")
        prev_rs = rs

    print()
    print("Note: '<<<' marks a sharp transition (bifurcation point).")
    print()


if __name__ == "__main__":
    main()
