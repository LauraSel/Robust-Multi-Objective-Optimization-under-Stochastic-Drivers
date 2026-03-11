# Robust-Multi-Objective-Optimization-under-Stochastic-Drivers
Robust Multi-Objective Optimization under Stochastic Drivers for Environmental Decision Support
This repository contains the MATLAB code used to generate the results reported in the paper:

**Robust Multi-Objective Optimization under Stochastic Drivers for Environmental Decision Support**

## Overview

The repository implements a simulation-optimization workflow for uncertainty-aware environmental decision support. The workflow combines:

- a deterministic stock-and-flow model (SFM),
- nominal multi-objective optimization,
- analysis of optimizer-induced variability across random seeds,
- prior-predictive robustness analysis under stochastic exogenous drivers.

The implementation is organized into three main steps:

1. **Step A — Nominal multi-objective optimization**
   - repeated runs of `gamultiobj`
   - construction of per-seed nondominated sets
   - candidate-set consolidation
   - preliminary hypervolume and spread diagnostics
     
2. **Reference-point construction**
   - computation of a single fixed hypervolume reference point from the nominal results
   - export of the reference point to `REF_global.mat`

3. **Step B — Optimizer variability analysis**
   - hypervolume computation with a fixed reference point
   - empirical attainment analysis
   - construction of the final candidate set for robustness analysis

4. **Step C — Prior-predictive robustness analysis**
   - Monte Carlo propagation of stochastic exogenous drivers
   - computation of robustness metrics
   - filtering and ranking of candidate interventions

## Software requirements

- MATLAB R2024b
- Global Optimization Toolbox
- Parallel Computing Toolbox (used for Monte Carlo robustness analysis with `parfor`)

## Main files

- `STEP_A.m` — nominal multi-objective optimization
- `MAKE_REF.m` — construction of the fixed hypervolume reference point from nominal outputs
- `STEP_B.m` — optimizer variability and candidate-set construction
- `STEP_C.m` — prior-predictive robustness analysis
- `obj_fun.m` — objective-function evaluator
- `sample_theta.m` — sampling of stochastic exogenous drivers from elicited prior distributions
- `default_theta.m` — baseline parameter configuration
- `hv3d_min.m` — hypervolume computation in minimization form

## Reproducibility

The workflow was executed in MATLAB R2024b on a machine equipped with an Intel Core i7 octa-core processor and 16 GB RAM.

To ensure reproducibility:

- repeated nominal runs use fixed random seeds,
- Monte Carlo simulations use fixed pseudo-random initialization,
- intermediate outputs are saved in `.mat` and `.csv` formats.

## How to run

A typical execution order is:

1. Run `STEP_A_nominal_MOO.m`
2. Run `MAKE_REF.m`
3. Run `STEP_B_variability.m`
4. Run `STEP_C_robustness.m`

## Data availability

This repository contains code only. Any case-study-specific inputs that cannot be redistributed should be described separately here.  
If all inputs are shareable, add them to the repository and document them in this section.

## Citation

If you use this code, please cite the associated paper.

## License

This repository is released under the MIT License.
