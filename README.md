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

---

## Supplementary archive: `validation_pipelines.zip`

Independent **transcriptomic validation** of the BOLD framework, corresponding
to Results §3.4 and Table 1 of the manuscript. This archive is self-contained
and complements the dynamical-modelling code above: it tests the model's
pre-specified qualitative predictions in three public GEO datasets, using the
same three-variable logic (X = inflammation, Y = osteogenesis, Z = metabolic
fitness) and the identical composite score formula
`RegScore = 0.4·Z + 0.4·Y + 0.2·(1 − X)`.

### Datasets and predictions tested

| Dataset | Assay | Prediction tested |
|---|---|---|
| **GSE234451** | Murine fracture-healing snRNA-seq (13 502 nuclei, Day 0/3/5/7) | Metabolic-axis attractor architecture; inflammatory-axis architecture; temporal escape trajectory; cell-type partitioning |
| **GSE240390** | Murine fracture-callus bulk RNA-seq (30 samples, Lean vs DIO, Day 3–21) | Pathological checkpoint impairment (Day-7 inflammatory-to-osteogenic transition) |
| **GSE180504** | Human BMSC bulk RNA-seq (Lean / Obese / T2D, n = 2 per group) | Pathological attractor deepening (T2D lowers the regenerative score) |

### Contents

| Path | Purpose |
|---|---|
| `README.md` | Full documentation: dataset table, quick-start, output inventory, reproduced-results table, environment, implementation notes |
| `run_all.sh` | Master runner — executes pipelines 01→04 in order, tees stdout to `logs/`; supports `DATA_DIR` / `OUTDIR` overrides |
| `config/gene_modules.py` | Single source of truth for all gene sets: X/Y/Z modules (mouse + human orthologues), M1/M2 markers, cell-type markers, Ensembl→symbol map, sample↔time/condition maps |
| `scripts/common.py` | Shared helpers: `minmax()`, `regscore()`, `module_score()`, `cohens_d()`, deterministic JSON writer |
| `scripts/01_gse234451_snrna.py` | GSE234451 snRNA-seq pipeline: time-point annotation, composite score, M1-like/M2-like classification, cell-type partitioning, marker-gene trajectories (P1/P2/P3/P7) |
| `scripts/02_gse240390_bulk.py` | GSE240390 bulk RNA-seq pipeline: Ensembl→symbol mapping, log2(CPM+1) module scoring, per-day Lean-vs-DIO statistics, pre-specified Day-7 checkpoint test (P5) |
| `scripts/03_gse180504_human_bmsc.py` | GSE180504 human BMSC pipeline: human-orthologue module scoring, group-level RegScore, directional pairwise comparisons (P4) |
| `scripts/04_validation_figure.py` | Assembles all pipeline outputs into the 6-panel validation figure |
| `outputs/` | All results committed: per-dataset CSVs, machine-readable `validation_summary.json` files, and `Figure_validation_pipeline.{png,pdf}` |
| `logs/` | Per-pipeline stdout logs providing the provenance trail for every number in Table 1 |

### Reproduced validation results

| BOLD prediction | Dataset | Result | Agreement |
|---|---|---|---|
| Metabolic-axis attractor architecture | GSE234451 | M2-like Z = +0.474 vs M1-like −0.121 (Δ = +0.595, d = 3.32) | ✓ |
| Inflammatory-axis architecture | GSE234451 | X did **not** segregate as predicted (Δ = +0.108, d = +0.56) | ✗ (reported explicitly) |
| Temporal escape trajectory | GSE234451 | X peaks at Day 5 (+0.064) and resolves by Day 7 (−0.013); RegScore 0.286 → 0.340 | ✓ |
| Cell-type partitioning | GSE234451 | Macrophages highest X; osteoblasts highest Y (Day 7 = 0.540) | ✓ |
| Pathological checkpoint impairment | GSE240390 | Day-7 ΔRegScore = −0.118; one-sided p = 0.050 (two-sided 0.101) | ✓ (preliminary) |
| Pathological attractor deepening | GSE180504 | T2D 0.341 vs Lean 0.472; Obese 0.551 (non-predicted direction) | ✓ (directional, n = 2) |

### Usage

```bash
unzip validation_pipelines.zip
cd validation_pipelines
bash run_all.sh
```

Requires `numpy`, `pandas`, `scipy`, `scanpy` (pipeline 01 only) and
`matplotlib` (figure step). All outputs are fully deterministic given the
input files; no random seeds are used. The negative inflammatory-axis result
is reproduced faithfully rather than tuned to force agreement.
