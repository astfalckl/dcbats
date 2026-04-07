# AR(2) Simulation and DC-BATS Reproducibility

This document provides a minimal, end-to-end set of shell commands to reproduce the AR(2) simulation, posterior estimation, Wasserstein aggregation, and final figure.

---

## Requirements

- R (≥ 4.2 recommended)
- R packages:
  - cmdstanr
  - posterior
  - optparse
  - ggplot2
  - dplyr
  - tidyr
- CmdStan installed and configured via cmdstanr

---

## Directory setup

mkdir -p figure_data \
mkdir -p figure_results \
mkdir -p case_study_pilot

---

## Notes

- For quick tests, reduce:
  - --T
  - --iter-warmup
  - --iter-sampling

# Make Figure 1

This provides the workflow of a single end-to-end simulation that reproduces Figure~1. Run all commands in terminal from the current directory.

---

## Step 1 — Simulate data

```bash
Rscript simulate_data.R \
  --T 10000 \
  --p 50 \
  --phi1 0.4 \
  --phi2 -0.6 \
  --output figure_data/ar2_sim.rds
```

---

## Step 2 — Fit full and subset posteriors

``` bash
Rscript fit_posteriors.R \
  --data figure_data/ar2_sim.rds \
  --stan-file model_ar2_errors.stan \
  --output-dir figure_results \
  --K 10,20 \
  --chains 4 \
  --parallel-chains 4 \
  --iter-warmup 1000 \
  --iter-sampling 1000
```

Outputs:

- figure_results/full_posterior_draws.rds
- figure_results/subset_posteriors_K10.rds
- figure_results/subset_posteriors_K20.rds

---

## Step 3 — Wasserstein aggregation

``` bash
Rscript wasserstein_average.R \
  --input figure_results/subset_posteriors_K10.rds \
  --output figure_results/wasserstein_beta_K10.rds

Rscript wasserstein_average.R \
  --input figure_results/subset_posteriors_K20.rds \
  --output figure_results/wasserstein_beta_K20.rds
```

---

## Step 4 — Generate figure

``` bash
Rscript make_figure_paper.R \
  --data figure_data/ar2_sim.rds \
  --full figure_results/full_posterior_draws.rds \
  --w10 figure_results/wasserstein_beta_K10.rds \
  --w20 figure_results/wasserstein_beta_K20.rds \
  --n-show 10 \
  --width 8 \
  --height 3 \
  --output figure_results/CI_ar2_errors_main.pdf
```

---

## Output

figure_results/CI_ar2_errors_main.pdf

---

# Frequentist Coverage (Table 1)

Please note that this takes a very long time. For minimal reproduction set
  - N_REPS=3
  - T_VAL=10000

in run_case_study_batch.sh

## Batch run simulations 

```bash
chmod +x run_case_study_batch.sh

./run_case_study_batch.sh
```

## Aggregate results into table

``` bash
Rscript summarise_coverage_case_study.r \
  --root case_study_pilot \
  --n-reps 50
```
