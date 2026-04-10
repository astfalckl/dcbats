# DC-BATS Code and Reproducibility

This repository contains the three case studies from the paper "Scalable Bayesian Inference for Time Series via Divide and Conquer". The folders are intended to be self-contained analysis directories rather than a single packaged R project.

## Structure

- `linear_regression/`: AR(2) linear regression simulation, Figure 1 workflow, and the repeated coverage study for Table 1
- `artfima/`: ARTFIMA simulation study, Wasserstein summaries, and example plots
- `air_quality/`: LA air quality case study for Table 3

## Minimal setup

- R
- folder-specific R package requirements are listed in each subdirectory README
- `cmdstanr` with a working CmdStan installation is needed for `linear_regression/` and `air_quality/`

## Notes

- The documentation is intentionally lightweight: enough to rerun the code, not a fully packaged workflow.
- Some simulations were originally run at larger scale on HPC. Where that matters, the subdirectory README points to smaller settings for a quick reproducibility check.
