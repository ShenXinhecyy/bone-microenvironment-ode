#!/usr/bin/env python3
"""
Monte Carlo robustness analysis of the bifurcation threshold.

Each of the 20 model parameters is independently perturbed by +/-20%
(uniform). For each of 10,000 perturbed parameter sets the saddle-node fold
k_ROS* is located by equilibrium continuation (Newton solve + bisection).

A parameter set is counted as "maintaining the qualitative prediction" if
(i) a high-inflammation attractor exists at k_ROS = 0 and (ii) a fold exists
within k_ROS in (0, 0.9).

Results (random seed 42, n = 10,000):
    - qualitative prediction maintained in 100.0% of sets
    - mean k_ROS* = 0.333, SD = 0.125
    - 95% percentile interval: 0.12-0.59
    - fold location dominated by alpha_X (Spearman rho = +0.64) and
      K_X (rho = -0.65)

Outputs: monte_carlo_folds.csv (one fold per valid parameter set).
Runtime: ~2-3 minutes on a desktop CPU.
"""

import numpy as np
import pandas as pd
from scipy import stats
from hysteresis_analysis import fold_forward, steady
from ode_model import DEFAULT_PARAMS_DM


def main(n_iter=10000, seed=42, out='monte_carlo_folds.csv'):
    rng = np.random.default_rng(seed)
    keys = list(DEFAULT_PARAMS_DM.keys())
    U = rng.uniform(0.8, 1.2, size=(n_iter, len(keys)))

    folds = np.full(n_iter, np.nan)
    for i in range(n_iter):
        params = {k: DEFAULT_PARAMS_DM[k] * U[i, j] for j, k in enumerate(keys)}
        try:
            folds[i] = fold_forward(params)
        except Exception:
            pass

    valid = ~np.isnan(folds)
    f = folds[valid]
    print(f"Qualitative prediction maintained: {valid.mean() * 100:.1f}% "
          f"({valid.sum()}/{n_iter})")
    print(f"mean k_ROS* = {f.mean():.3f}, SD = {f.std():.3f}")
    print(f"95% percentile interval: {np.percentile(f, 2.5):.2f}-"
          f"{np.percentile(f, 97.5):.2f}")
    W, p = stats.shapiro(f[:5000])
    print(f"Shapiro-Wilk W = {W:.3f}, p = {p:.4f}")

    print("\nParameter influence on fold location (Spearman rank correlation):")
    sens = sorted(
        ((k, stats.spearmanr(U[valid, j], f).statistic)
         for j, k in enumerate(keys)),
        key=lambda x: -abs(x[1]))
    for k, rho in sens[:6]:
        print(f"  {k:16s} rho = {rho:+.3f}")

    pd.DataFrame({'k_ROS_star': f}).to_csv(out, index=False)
    print(f"\nSaved {len(f)} fold values to {out}")


if __name__ == '__main__':
    main()
