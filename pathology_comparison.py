#!/usr/bin/env python3
"""
Cross-pathology comparison: simulate Healthy, Diabetes, Osteoporosis,
and Osteomyelitis conditions, with and without therapeutic intervention.
"""

import numpy as np
from ode_model import DEFAULT_PARAMS_DM, simulate, get_final_state, regscore


def main():
    # Build pathology-specific parameter sets
    params_healthy = DEFAULT_PARAMS_DM.copy()
    params_healthy['alpha_X'] = 0.25  # Low inflammation drive

    params_dm = DEFAULT_PARAMS_DM.copy()  # Diabetes: high alpha_X

    params_op = DEFAULT_PARAMS_DM.copy()
    params_op['alpha_X'] = 0.55
    params_op['alpha_Z'] = 0.03  # Impaired angiogenesis

    params_om = DEFAULT_PARAMS_DM.copy()
    params_om['alpha_X'] = 0.90  # Extreme inflammation
    params_om['gamma_ZX_damage'] = 0.40  # Severe damage coupling

    pathologies = [
        ('Healthy', params_healthy, [0.2, 0.6, 0.7]),
        ('Diabetes', params_dm, [0.7, 0.15, 0.20]),
        ('Osteoporosis', params_op, [0.7, 0.15, 0.20]),
        ('Osteomyelitis', params_om, [0.7, 0.15, 0.20]),
    ]

    # Therapeutic intervention (Triple Push-Pull)
    therapy = {'k_ROS': 0.5, 'alpha_AMPK': 0.35, 'E_matrix': 0.2}
    no_therapy = {'k_ROS': 0.0, 'alpha_AMPK': 0.0, 'E_matrix': 0.0}

    print("=" * 100)
    print("Cross-Pathology Comparison")
    print("=" * 100)
    print(f"{'Pathology':15s} {'Before (X,Y,Z,RS)':>40s} {'After (X,Y,Z,RS)':>40s} {'ΔRS':>8s}")
    print("-" * 105)

    for pname, pparams, init in pathologies:
        sol_before = simulate(pparams, no_therapy, init)
        sol_after = simulate(pparams, therapy, init)

        Xb, Yb, Zb = get_final_state(sol_before)
        Xa, Ya, Za = get_final_state(sol_after)
        rs_b = regscore(Xb, Yb, Zb)
        rs_a = regscore(Xa, Ya, Za)

        print(f"{pname:15s}  X={Xb:.3f} Y={Yb:.3f} Z={Zb:.3f} RS={rs_b:.3f}   "
              f"X={Xa:.3f} Y={Ya:.3f} Z={Za:.3f} RS={rs_a:.3f}   {rs_a - rs_b:+7.3f}")

    print()


if __name__ == "__main__":
    main()
