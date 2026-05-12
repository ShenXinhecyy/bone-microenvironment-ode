#!/usr/bin/env python3
"""
Compare different therapeutic strategies for the Diabetes Mellitus model.
"""

import numpy as np
from ode_model import DEFAULT_PARAMS_DM, simulate, get_final_state, regscore


def main():
    params = DEFAULT_PARAMS_DM
    init_state = [0.7, 0.15, 0.20]

    strategies = [
        ('None (control)', {'k_ROS': 0.0, 'alpha_AMPK': 0.0, 'E_matrix': 0.0}),
        ('ROS clearance only', {'k_ROS': 0.5, 'alpha_AMPK': 0.0, 'E_matrix': 0.0}),
        ('AMPK activation only', {'k_ROS': 0.0, 'alpha_AMPK': 0.35, 'E_matrix': 0.0}),
        ('Stiffness only', {'k_ROS': 0.0, 'alpha_AMPK': 0.0, 'E_matrix': 0.2}),
        ('ROS + AMPK', {'k_ROS': 0.5, 'alpha_AMPK': 0.35, 'E_matrix': 0.0}),
        ('Triple Push-Pull', {'k_ROS': 0.5, 'alpha_AMPK': 0.35, 'E_matrix': 0.2}),
    ]

    print("=" * 70)
    print("Therapeutic Strategy Comparison (Diabetes Model)")
    print("=" * 70)
    print(f"{'Strategy':25s}  {'X':>8s}  {'Y':>8s}  {'Z':>8s}  {'RS':>8s}")
    print("-" * 60)

    for sname, control in strategies:
        sol = simulate(params, control, init_state)
        X, Y, Z = get_final_state(sol)
        rs = regscore(X, Y, Z)
        print(f"  {sname:25s}  {X:>7.4f}  {Y:>7.4f}  {Z:>7.4f}  {rs:>7.4f}")

    print()


if __name__ == "__main__":
    main()
