#!/usr/bin/env python3
"""
Hysteresis / permanence analysis of the BOLD model (diabetes parameterisation).

Forward fold (disappearance of the inflammatory attractor as k_ROS increases):
    k_ROS* ~ 0.333 (equilibrium continuation).

Reverse fold (disappearance of the regenerative attractor as k_ROS decreases):
    no reverse fold exists for any k_ROS >= -0.3, i.e. the regenerative
    attractor persists even at k_ROS = 0 (and below zero, an unphysical
    ROS-generating regime). Interpretation: once the system has crossed the
    forward fold, withdrawal of the intervention does NOT return it to the
    inflammatory attractor - the regenerative state transition is permanent
    within the model. The diabetic baseline itself is bistable; the pathology
    is a problem of basin occupancy, not of attractor existence.
"""

import numpy as np
from bifurcation_scan import steady
from ode_model import DEFAULT_PARAMS_DM, regscore


def fold_forward(params=DEFAULT_PARAMS_DM):
    def infl(k):
        s = steady({'k_ROS': k}, [0.79, 0.14, 0.06], params)
        return s is not None and s[0] > 0.4
    lo, hi = 0.0, 0.9
    if not infl(lo):
        return np.nan
    for _ in range(45):
        mid = (lo + hi) / 2
        if infl(mid):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def fold_reverse(params=DEFAULT_PARAMS_DM, k_min=-0.5):
    """Returns the reverse fold, or None if the regenerative attractor
    persists down to k_min (i.e. no physiological reverse fold)."""
    def regen(k):
        s = steady({'k_ROS': k}, [0.05, 0.80, 0.95], params)
        return s is not None and s[2] > 0.8
    if regen(k_min):
        return None
    lo, hi = k_min, 0.9
    for _ in range(45):
        mid = (lo + hi) / 2
        if regen(mid):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


if __name__ == '__main__':
    ff = fold_forward()
    fr = fold_reverse()
    print(f"Forward fold (inflammatory branch annihilation): k_ROS* = {ff:.4f}")
    if fr is None:
        print("Reverse fold: none found down to k_ROS = -0.5.")
        print("=> The regenerative attractor persists at k_ROS = 0; above-fold")
        print("   interventions produce a PERMANENT state transition (no relapse")
        print("   upon intervention withdrawal within the model horizon).")
        s0 = steady({'k_ROS': 0.0}, [0.05, 0.80, 0.95])
        print(f"   Regenerative attractor at k_ROS=0: X={s0[0]:.3f} "
              f"Y={s0[1]:.3f} Z={s0[2]:.3f} RS={regscore(*s0):.3f}")
    else:
        print(f"Reverse fold: k_ROS = {fr:.4f} (hysteresis width "
              f"{ff - fr:.4f})")
