# ARTFIMA Simulation

## Overview

This repository implements a fully reproducible pipeline for:

1.  Simulating ARTFIMA time series
2.  Fitting posterior distributions via a custom Metropolis--Hastings
    sampler
3.  Constructing Wasserstein barycenter approximations using subset
    posteriors
4.  Aggregating results across Monte Carlo replicates

------------------------------------------------------------------------

## Some notes

Some code is borrowed from an archived version of artfima. The package
is no longer available on CRAN, but is copied in this repo, and the required
functions are hard-called (i.e. no installation necessary).

It is also a fairly lengthy amount of simulation. For reproducibility I would
set 

------------------------------------------------------------------------

## Directory Structure

    project/
      R/
        helpers_artfima.R
        fast_acf_from_psd.R
        mh_artfima.R
        wasserstein_utils.R

      artfima/  <--- [archived version of artfima]

      scripts/
        simulate_data.R
        fit_posteriors.R
        wasserstein_average.R
        aggregate_results.R

      data/
        simulated/

      results/
        fits/
        summaries/

------------------------------------------------------------------------

## Pipeline Description

### 1. Simulation (`simulate_data.R`)

-   Generates ARTFIMA time series using a fast PSD → ACF calculation
-   Outputs:
    -   data/simulated/<scenario>/<sim_id>/data.rds
    -   data/simulated/<scenario>/<sim_id>/meta.rds

------------------------------------------------------------------------

### 2. Posterior Fitting (`fit_posteriors.R`)

-   Fits:
    -   Full posterior (temper = 1)
    -   Subset posteriors for K = 5, 10, 20 (tempered likelihood)
-   Outputs:
    -   results/fits/<scenario>/<sim_id>/full.rds
    -   results/fits/<scenario>/<sim_id>/K05/k01.rds

------------------------------------------------------------------------

### 3. Wasserstein Averaging (`wasserstein_average.R`)

-   Converts posterior draws into empirical quantile functions: Q(u), u
    ∈ {0.01,...,0.99}
-   Computes subset posterior barycenters: average of quantiles across
    subsets
-   Outputs:
    -   results/summaries/posterior_quantiles.rds

------------------------------------------------------------------------

### 4. Aggregation (`aggregate_results.R`)

-   Computes relative Wasserstein distances
-   Aggregates across simulation replicates
-   Outputs:
    -   results/summaries/wasserstein_by_sim.rds
    -   results/summaries/wasserstein_summary.rds

------------------------------------------------------------------------

## Execution

Run the pipeline sequentially:

``` bash
Rscript scripts/simulate_data.R
Rscript scripts/fit_posteriors.R
Rscript scripts/wasserstein_average.R
Rscript scripts/aggregate_results.R
```

------------------------------------------------------------------------

## Session Info

```r
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.4

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Australia/Sydney
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tibble_3.3.1     readr_2.2.0      purrr_1.2.1      stringr_1.6.0    tidyr_1.3.2      dplyr_1.2.0      SuperGauss_2.0.4
 [8] gsl_2.1-9        devtools_2.5.0   usethis_3.2.1    ltsa_1.4.6.1    

loaded via a namespace (and not attached):
 [1] jsonlite_2.0.0    compiler_4.5.2    tidyselect_1.2.1  Rcpp_1.1.1        fastmap_1.2.0     R6_2.6.1          generics_0.1.4   
 [8] pillar_1.11.1     tzdb_0.5.0        rlang_1.1.7       utf8_1.2.6        cachem_1.1.0      stringi_1.8.7     fs_2.0.1         
[15] fftw_1.0-9        pkgload_1.5.0     memoise_2.0.1     cli_3.6.5         withr_3.0.2       magrittr_2.0.4    hms_1.1.4        
[22] lifecycle_1.0.5   vctrs_0.7.1       glue_1.8.0        sessioninfo_1.2.3 pkgbuild_1.4.8    tools_4.5.2       pkgconfig_2.0.3  
[29] ellipsis_0.3.2 
```