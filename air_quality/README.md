# LA Air Quality

## Overview

This repository implements a fully reproducible pipeline for the LA particulate matter case study given in Section 5. As runs are not iterated over repeated simulation as in the simulations in Section 4, the pipeline is much simpler.

    project/
        dcc_garch.stan
        la_pm_log_positive.csv
        run_analysis.r

This runs in reasonable time (~1 hr) on a standard laptop, but if you want to expedite the results, cut down `iter_warmup` and `iter_warmup` to 1000.

The code used `cmdstanr` which is not available on CRAN and must be installed as

```r
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```

------------------------------------------------------------------------

## Running the analysis

The full analysis is contained in `run_analysis.r`. It is simple enough and can be housed in a single script. Table~3 is exported as `results_intervals.csv`.

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