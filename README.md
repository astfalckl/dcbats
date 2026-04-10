# DC-BATS Code and Reproducibility

This repository contains the three case studies from the paper "Scalable Bayesian Inference for Time Series via Divide and Conquer". The folders are intended to be self-contained analysis directories rather than a single packaged R project.

## Structure

- `linear_regression/`: AR(2) linear regression simulation, Figure 1 workflow, and the repeated coverage study for Table 1
- `artfima/`: ARTFIMA simulation study, Wasserstein summaries, and the Figure 2 example plot
- `air_quality/`: LA air quality case study for Table 3

## Results map

- Figure 1: `linear_regression/` via `make_figure_paper.r`; main output `figure_results/CI_ar2_errors_main.pdf`
- Table 1: `linear_regression/` via `run_case_study_batch.sh` and `summarise_coverage_case_study.r`; summary output `case_study/coverage_summary_table.csv`
- Table 2: `artfima/` via `scripts/aggregate_results.r`; summary outputs `results/summaries/wasserstein_summary.rds` and `results/summaries/wasserstein_summary_long.rds`
- Figure 2: `artfima/` via `scripts/plot_examples.r`; set `sim_id_target <- "012"` to match the paper figure
- Table 3: `air_quality/` via `run_analysis.r`; output `results_intervals.csv`

## Minimal setup

- R
- folder-specific R package requirements are listed in each subdirectory README
- `cmdstanr` with a working CmdStan installation is needed for `linear_regression/` and `air_quality/`

## Notes

- The documentation is intentionally lightweight: enough to rerun the code, not a fully packaged workflow.
- Some simulations were originally run at larger scale on HPC. Where that matters, the subdirectory README points to smaller settings for a quick reproducibility check.
- The repository documents the checked-in scripts used to generate the main results; cluster-specific submission files are not included.
