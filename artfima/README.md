# ARTFIMA Simulation

This folder contains the ARTFIMA simulation study and the Figure 2 example plot used in the paper.

## Results in paper

- Table 2: `scripts/aggregate_results.r`
- Figure 2: `scripts/plot_examples.r`

## Requirements

- R
- `ltsa`, `devtools`, `gsl`, `SuperGauss`, `dplyr`, `tidyr`, `stringr`, `purrr`, `readr`, `tibble`

An archived copy of the `artfima` code is vendored in the nested `artfima/` directory, so there is no separate installation step for that dependency.

## Structure

- `R/`: helper functions, MH sampler code, and Wasserstein utilities
- `scripts/simulate_data.r`: simulate datasets under `data/simulated/`
- `scripts/fit_posteriors.r`: fit full and subset posteriors under `results/fits/`
- `scripts/wasserstein_average.r`: compute posterior quantile summaries
- `scripts/aggregate_results.r`: aggregate relative Wasserstein distances across simulations
- `scripts/plot_examples.r`: generate an example density comparison plot under `results/summaries/`

## Run

Run the pipeline from this directory:

```bash
Rscript scripts/simulate_data.r
Rscript scripts/fit_posteriors.r
Rscript scripts/wasserstein_average.r
Rscript scripts/aggregate_results.r
Rscript scripts/plot_examples.r
```

The summary outputs are written under `results/summaries/`.

## Notes

- The full simulation study uses `n_simulations <- 100` in `scripts/simulate_data.r` and takes a long time.
- For a lighter reproducibility run, set `n_simulations <- 2` before running the pipeline.
- Not all generated simulation outputs are checked into the repository because of size.
- The paper-scale version of this study was run in parallel on HPC, but the repository only includes the analysis scripts themselves, not cluster-specific `.pbs` submission files.
- `scripts/plot_examples.r` uses a fixed simulation index. Set `sim_id_target <- "012"` to match the paper figure noted in the script comments; the checked-in default can be changed if you want a different replicate.
