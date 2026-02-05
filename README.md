# Copula-Based Joint Regression with Beta Margins: Simulation Code

This repository contains reproducible simulation code for "Joint Copula-Based Beta Regression with
Covariate-Dependent Dispersion for Malaria
Indicators". The code supports the simulation studies reported in the accompanying manuscript and is designed to evaluate finite-sample performance, numerical stability, and convergence behavior of likelihood-based joint models for bounded outcomes.

## Overview

The simulation study is motivated by applications involving bounded continuous responses and focuses on assessing model performance under varying data-generating scenarios. In particular, the simulations examine:

- Different sample sizes ($n = 400, 800, 1600$)
- Varying levels of cross-outcome dependence
- Nonlinear covariate effects in the marginal mean structures
- Constant versus covariate-dependent dispersion specifications
- Numerical stability and convergence behavior of the estimation procedure

The joint dependence structure is modeled using a copula, while marginal distributions are specified using Beta regression with flexible mean and dispersion components.

## Repository Structure

- `R/`  
  Contains R scripts for data generation, model fitting, and extraction of simulation results.

- `simulation/`  
  Scripts for running simulation scenarios across different sample sizes, dependence levels, and dispersion structures.

- `results/`  
  Output files summarizing parameter estimates, convergence rates, and performance metrics (e.g., bias, RMSE).

- `utils/`  
  Helper functions used throughout the simulation study.

## Software Requirements

The simulations are implemented in **R** (version ≥ 4.1.0) and rely on the following packages:

- `mgcv` — penalized regression splines and REML-based smoothing parameter selection  
- `gamlss` — distributional regression components  
- `copula` — copula construction and likelihood evaluation  
- `GJRM` — joint regression modeling framework  
- `tidyverse` — data manipulation and visualization

All packages are available from CRAN.

## Reproducibility

Simulation scripts are fully reproducible. Each script sets a random seed and documents the data-generating mechanisms used in the study. The results reported in the manuscript were obtained by running the simulation scripts provided in this repository without modification.

To reproduce the simulation results:

1. Install the required R packages.
2. Run the main simulation scripts in the `simulation/` directory.
3. Aggregate and summarize results using the scripts in the `results/` directory.

## Notes on Computation

All simulations were conducted using likelihood-based estimation. For the sample sizes considered, model fitting was computationally feasible and did not pose practical difficulties. Occasional non-convergence may occur in smaller samples when dependence parameters are weakly identified; such cases are documented in the simulation output and excluded from performance summaries, consistent with standard practice.

## Citation

If you use this code or build upon it in your own work, please cite the associated manuscript.

## Contact

For questions or comments regarding the code, please contact:



