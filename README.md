# Bone Microenvironment ODE Model

A three-variable ODE model describing the bone regeneration microenvironment, featuring X (inflammation/ROS), Y (osteogenesis), and Z (angiogenesis) dynamics with bistability and bifurcation behavior.

## Model Equations

The model uses a unified saturation mechanism (Model C — the final corrected version):

```
dX/dt = αX·Hill(X)·(1−X) − βX·X − γXY·Hill(Y)·X − γZX·Hill(Z)·X − kROS·X
dY/dt = (αYb + γYZ·Z + βY·Y + Emat)·(1−Y) − κXY·Hill(X)·Y − δY·Y
dZ/dt = −γZXd·X·Z + (αZ·(1−X) + βZ·Z + αAMPK)·(1−Z)
```

### Variables
- **X** — Inflammation/ROS level (0–1)
- **Y** — Osteogenesis/bone formation activity (0–1)
- **Z** — Angiogenesis/vascularization (0–1)

### Control Parameters
- `k_ROS` — ROS clearance rate (therapeutic)
- `alpha_AMPK` — AMPK activation (metabolic rescue)
- `E_matrix` — Extracellular matrix stiffness (mechanical support)

### RegScore (Regeneration Score)
```
RegScore = 0.4·Z + 0.4·Y + 0.2·(1−X)
```

## Files

| File | Description |
|------|-------------|
| `ode_model.py` | Core ODE system with Hill function and unified saturation |
| `bistability_test.py` | Test for bistability across initial conditions |
| `bifurcation_scan.py` | Bifurcation analysis with k_ROS as control parameter |
| `pathology_comparison.py` | Cross-pathology comparison (Healthy, Diabetes, Osteoporosis, Osteomyelitis) |
| `strategy_comparison.py` | Therapeutic strategy comparison (ROS, AMPK, Stiffness) |
| `requirements.txt` | Python dependencies |

## Usage

```bash
pip install -r requirements.txt
python bistability_test.py
python bifurcation_scan.py
python pathology_comparison.py
python strategy_comparison.py
```

## Model Versions

Three versions were developed and compared:
- **Model A** — Original (no saturation terms, unbounded risk)
- **Model B** — Paper version (partial saturation)
- **Model C** — Unified saturation (final, most stable)

Model C is the recommended version.
