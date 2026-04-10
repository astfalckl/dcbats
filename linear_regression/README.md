# Linear Regression with AR(2) Errors

This folder contains the AR(2) linear regression example used for Figure 1 and the repeated simulation used for Table 1.

## Requirements

- R
- `cmdstanr`, `posterior`, `optparse`, `ggplot2`, `dplyr`, `tidyr`
- CmdStan installed and configured via `cmdstanr`

## Reproduce Figure 1

Run the following from this directory:

```bash
mkdir -p figure_data figure_results

Rscript simulate_data.r \
  --T 10000 \
  --p 50 \
  --phi1 0.4 \
  --phi2 -0.6 \
  --output figure_data/ar2_sim.rds

Rscript fit_posteriors.r \
  --data figure_data/ar2_sim.rds \
  --stan-file model_ar2_errors.stan \
  --output-dir figure_results \
  --K 10,20 \
  --chains 4 \
  --parallel-chains 4 \
  --iter-warmup 1000 \
  --iter-sampling 1000

Rscript wasserstein_average.r \
  --input figure_results/subset_posteriors_K10.rds \
  --output figure_results/wasserstein_beta_K10.rds

Rscript wasserstein_average.r \
  --input figure_results/subset_posteriors_K20.rds \
  --output figure_results/wasserstein_beta_K20.rds

Rscript make_figure_paper.r \
  --data figure_data/ar2_sim.rds \
  --full figure_results/full_posterior_draws.rds \
  --w10 figure_results/wasserstein_beta_K10.rds \
  --w20 figure_results/wasserstein_beta_K20.rds \
  --n-show 10 \
  --width 8 \
  --height 3 \
  --output figure_results/CI_ar2_errors_main.pdf
```

The main output is `figure_results/CI_ar2_errors_main.pdf`.

## Coverage study

`run_case_study_batch.sh` runs the repeated simulation study and writes outputs under `case_study/`.

The paper-scale settings in that script (`N_REPS=50`, `T_VAL=100000`) are slow and were intended for HPC use. The repository includes the batch driver script itself, but not scheduler-specific `.pbs` submission files. For a lighter reproducibility run, reduce those values first, then run:

```bash
chmod +x run_case_study_batch.sh
./run_case_study_batch.sh

Rscript summarise_coverage_case_study.r \
  --root case_study \
  --n-reps {N_REPS}
```

This writes `case_study/coverage_summary_long.csv` and `case_study/coverage_summary_table.csv`.

## Notes

- For quick tests, the main parameters to reduce are `--T`, `--iter-warmup`, and `--iter-sampling`.
- If you want to rerun the larger study on a cluster, `run_case_study_batch.sh` is the entry point to wrap in your local scheduler.
