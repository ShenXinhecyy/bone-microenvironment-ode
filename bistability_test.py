#!/usr/bin/env python3
"""
Test for bistability: run the model from different initial conditions
and check if multiple distinct steady states exist.
"""

import numpy as np
from ode_model import DEFAULT_PARAMS_DM, simulate, get_final_state, regscore


def main():
    params = DEFAULT_PARAMS_DM
    control = {'k_ROS': 0.0, 'alpha_AMPK': 0.0, 'E_matrix': 0.0}

    initial_conditions = [
        [0.7, 0.15, 0.2],   # High inflammation start
        [0.5, 0.3, 0.4],    # Moderate start
        [0.3, 0.5, 0.6],    # Low inflammation start
        [0.1, 0.8, 0.8],    # Healing start
    ]

    print("=" * 70)
    print("Bistability Test — Bone Microenvironment ODE Model")
    print("=" * 70)

    results = {}
    for init in initial_conditions:
        sol = simulate(params, control, init)
        X, Y, Z = get_final_state(sol)
        rs = regscore(X, Y, Z)
        bounds_ok = all(-0.01 <= v <= 1.01 for v in [X, Y, Z])
        status = "✓" if bounds_ok else "✗"
        results[tuple(init)] = (X, Y, Z, rs)
        print(f"  Init {init} → X={X:.4f}, Y={Y:.4f}, Z={Z:.4f}, RS={rs:.4f}  [{status}]")

    unique_rs = set(round(r[3], 2) for r in results.values())
    print(f"\n  Distinct steady states: {len(unique_rs)}")
    if len(unique_rs) > 1:
        print("  ✅ Bistability detected!")
    else:
        print("  ❌ No bistability (single attractor)")

    print()


if __name__ == "__main__":
    main()
