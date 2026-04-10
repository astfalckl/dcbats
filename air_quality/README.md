# LA Air Quality

This folder contains the LA particulate matter case study used for Table 3.

## Requirements

- R
- `cmdstanr`, `posterior`, `data.table`
- CmdStan installed and configured via `cmdstanr`

`cmdstanr` is not on CRAN. Install it with:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
```

## Files

- `dcc_garch.stan`: Stan model
- `la_pm_log_positive.csv`: input data
- `run_analysis.r`: end-to-end analysis script

## Run

From this directory, run:

```bash
Rscript run_analysis.r
```

This writes `results_intervals.csv`, which reproduces the interval summary reported in Table 3.

## Notes

- With the checked-in settings, the script fits the full model plus 10 subset models and is not especially quick.
- For a lighter reproducibility run, reduce `iter_warmup` and `iter_sampling` inside `run_analysis.r`.
