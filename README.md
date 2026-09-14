# BOLD Model Code — Bone Microenvironment Dynamical Model

Companion code for: "Breaking the Inflammatory Attractor: Quantitative Design
Rules for Smart Biomaterials from Integrated Evidence Synthesis and Dynamical
Modelling" (*Bone Research* submission).

The BOLD (Bistable Osteo-immune Logic for Design) model is a three-variable
nonlinear ODE system — X (pro-inflammatory activity), Y (pro-osteogenic
activity), Z (metabolic fitness) — with 20 parameters and 3 control inputs
(k_ROS, alpha_AMPK, E_matrix). Current parameterisation version: v5
(K_Z = 0.42; sole change from v4, placing the saddle-node fold at
k_ROS* = 0.333 by equilibrium continuation).

## Files

| File | Purpose |
|---|---|
| `ode_model.py` | Core ODE system, default diabetic parameter set (v5), simulator, RegScore |
| `bistability_test.py` | Confirms two distinct attractors from different initial states |
| `bifurcation_scan.py` | Three-branch equilibrium continuation (Newton + bisection); fold k_ROS* = 0.333; finite-time grid cross-check (0.338) |
| `hysteresis_analysis.py` | Forward/reverse continuation; shows no reverse fold down to k_ROS = -0.3 (permanence of the transition) |
| `monte_carlo_robustness.py` | 10 000-run ±20% uniform parameter perturbation (seed 42); fold maintained in 9 999/10 000 sets; Spearman rank sensitivity |
| `sensitivity_sweeps.py` | One-at-a-time control-authority sweeps (k_ROS 0.80 > alpha_AMPK 0.44 > E_matrix 0.18) |
| `strategy_comparison.py` | Push–Pull strategy comparison (none 0.12; ROS(0.5) 0.92; ROS(0.3)+AMPK 0.92 vs additive 0.72; Triple 0.94) |
| `pathology_comparison.py` | Cross-pathology gains: OM +0.84, DM +0.82, OP +0.80 |
| `monte_carlo_folds.csv` | Fold locations from the 9 999 valid Monte Carlo runs (Supplementary Data S2) |

## Quick start

```bash
pip install -r requirements.txt
python bistability_test.py
python bifurcation_scan.py
python hysteresis_analysis.py
python strategy_comparison.py
python pathology_comparison.py
python sensitivity_sweeps.py
python monte_carlo_robustness.py   # ~2-3 min
```

## Reproduced headline results (K_Z = 0.42)

- Bistability at diabetic baseline; inflammatory attractor RS = 0.12
- Saddle-node fold: equilibrium continuation 0.333; finite-time grid 0.338
- Reverse continuation: no fold above k_ROS = -0.3 → permanent transition
- Monte Carlo (n = 10 000, seed 42): mean k_ROS* = 0.333, SD = 0.125,
  95% percentile interval 0.12–0.59; dominant parameters K_X (ρ = -0.65)
  and alpha_X (ρ = +0.64)
- Sub-threshold synergy: ROS(0.3)+AMPK = 0.92 vs naive additive 0.72
