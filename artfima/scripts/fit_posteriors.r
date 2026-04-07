# ------------------------------------------------------------------------------
# 02_fit_posteriors.R
#
# Fit full-data and subset posteriors for previously simulated ARTFIMA datasets.
#
# This script:
#   1. loads shared helpers
#   2. reads simulated datasets from data/simulated/
#   3. fits:
#        - full posterior
#        - K = 5 subset posteriors (temper = 5)
#        - K = 10 subset posteriors (temper = 10)
#        - K = 20 subset posteriors (temper = 20)
#   4. saves chains + run metadata
#
# It does NOT compute Wasserstein summaries.
# ------------------------------------------------------------------------------

rm(list = ls())

# ------------------------------------------------------------------------------
# Load shared helpers
# ------------------------------------------------------------------------------

source("R/helpers_artfima.R")
load_required_packages()
source_artfima_code()

source("R/fast_acf_from_psd.R")
source("R/mh_artfima.R")

# ------------------------------------------------------------------------------
# User-configurable settings
# ------------------------------------------------------------------------------

simulation_root <- file.path("data", "simulated")
fit_root <- file.path("results", "fits")

ensure_project_subdirs(c(
  "results",
  fit_root
))

# MCMC controls
n_iter <- 5000
burn_in <- 1000

# Initial values:
# keep these aligned with the data-generating values, as in the current script
use_truth_as_init <- TRUE

# If FALSE, and fit files already exist, the script skips them
overwrite_existing <- FALSE

# Which subset partitions to fit
K_values <- c(5L, 10L, 20L)

# Whether to save data with the fit object metadata
store_data_with_fit <- TRUE

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

list_scenario_dirs <- function(root) {
  if (!dir.exists(root)) {
    stop(sprintf("Simulation root not found: %s", root), call. = FALSE)
  }

  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs[order(basename(dirs))]
}

list_sim_dirs <- function(scenario_dir) {
  dirs <- list.dirs(scenario_dir, recursive = FALSE, full.names = TRUE)
  dirs[order(basename(dirs))]
}

read_simulation_bundle <- function(sim_dir) {
  data_path <- file.path(sim_dir, "data.rds")
  meta_path <- file.path(sim_dir, "meta.rds")

  if (!file.exists(data_path)) {
    stop(sprintf("Missing simulation data file: %s", data_path), call. = FALSE)
  }

  if (!file.exists(meta_path)) {
    stop(sprintf("Missing simulation metadata file: %s", meta_path), call. = FALSE)
  }

  list(
    x = readRDS(data_path),
    meta = readRDS(meta_path)
  )
}

make_init_values <- function(meta, x) {
  if (isTRUE(use_truth_as_init)) {
    list(
      d = meta$scenario$d,
      lambda = meta$scenario$lambda,
      phi = meta$scenario$phi,
      sigma2 = 1
    )
  } else {
    list(
      d = 0.2,
      lambda = 0.01,
      phi = 0.5,
      sigma2 = stats::var(x)
    )
  }
}

fit_exists <- function(path) {
  file.exists(path)
}

save_fit_bundle <- function(
  fit_obj,
  out_path,
  run_info
) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  out <- list(
    draws = fit_obj,
    run_info = run_info
  )

  saveRDS(out, out_path)
  invisible(out_path)
}

fit_single_chain <- function(
  x,
  init,
  sd_prop,
  temper,
  n_iter,
  true_pars = NULL,
  subset_index = NULL,
  K = NULL,
  store_data = TRUE
) {
  fit_time <- system.time({
    draws <- mh_artfima(
      x = x,
      n_iter = n_iter,
      init = init,
      sd_prop = sd_prop,
      incl_sample = list(
        sample_d = TRUE,
        sample_lambda = FALSE,
        sample_phi = TRUE,
        sample_sigma2 = FALSE
      ),
      temper = temper
    )
  })

  run_info <- list(
    elapsed_seconds = unname(fit_time["elapsed"]),
    accept_rate = attr(draws, "accept_rate"),
    sd_prop = sd_prop,
    init = init,
    temper = temper,
    n_iter = n_iter,
    burn_in = burn_in,
    subset_index = subset_index,
    K = K,
    true_pars = true_pars
  )

  if (isTRUE(store_data)) {
    run_info$data <- x
  }

  list(
    draws = draws,
    run_info = run_info
  )
}

make_subset_indices <- function(n, K) {
  if (n %% K != 0) {
    stop(sprintf("Series length %d is not divisible by K = %d.", n, K), call. = FALSE)
  }

  m <- n / K

  lapply(seq_len(K), function(k) {
    start <- (k - 1L) * m + 1L
    end <- k * m
    start:end
  })
}

# ------------------------------------------------------------------------------
# Main fitting loop
# ------------------------------------------------------------------------------

scenario_dirs <- list_scenario_dirs(simulation_root)

if (length(scenario_dirs) == 0L) {
  stop("No scenario directories found under data/simulated/.", call. = FALSE)
}

for (scenario_dir in scenario_dirs) {
  scenario_id <- basename(scenario_dir)
  sim_dirs <- list_sim_dirs(scenario_dir)

  message("")
  message(sprintf("Scenario: %s", scenario_id))

  for (sim_dir in sim_dirs) {
    sim_id <- basename(sim_dir)

    message(sprintf("  Simulation: %s", sim_id))

    sim_bundle <- read_simulation_bundle(sim_dir)
    x <- sim_bundle$x
    meta <- sim_bundle$meta

    n <- length(x)

    true_pars <- meta$scenario
    init <- make_init_values(meta = meta, x = x)
    sd_prop <- default_mcmc_proposals(
      phi = true_pars$phi,
      d = true_pars$d
    )

    # --------------------------------------------------------------------------
    # Full posterior
    # --------------------------------------------------------------------------

    full_out_path <- file.path(
      fit_root,
      scenario_id,
      sim_id,
      "full.rds"
    )

    if (!overwrite_existing && fit_exists(full_out_path)) {
      message("    Full: exists, skipping")
    } else {
      message("    Full: fitting")

      fit_obj <- fit_single_chain(
        x = x,
        init = init,
        sd_prop = sd_prop,
        temper = 1,
        n_iter = n_iter,
        true_pars = true_pars,
        subset_index = NA_integer_,
        K = 1L,
        store_data = store_data_with_fit
      )

      fit_obj$run_info$scenario_id <- scenario_id
      fit_obj$run_info$sim_id <- sim_id
      fit_obj$run_info$source_data_path <- file.path(sim_dir, "data.rds")
      fit_obj$run_info$source_meta_path <- file.path(sim_dir, "meta.rds")
      fit_obj$run_info$fit_type <- "full"

      save_fit_bundle(
        fit_obj = fit_obj$draws,
        out_path = full_out_path,
        run_info = fit_obj$run_info
      )
    }

    # --------------------------------------------------------------------------
    # Subset posteriors
    # --------------------------------------------------------------------------

    for (K in K_values) {
      subset_dir <- file.path(
        fit_root,
        scenario_id,
        sim_id,
        sprintf("K%02d", K)
      )

      subset_indices <- make_subset_indices(n = n, K = K)

      message(sprintf("    K = %d", K))

      for (k in seq_len(K)) {
        subset_out_path <- file.path(
          subset_dir,
          sprintf("k%02d.rds", k)
        )

        if (!overwrite_existing && fit_exists(subset_out_path)) {
          message(sprintf("      subset %02d: exists, skipping", k))
          next
        }

        x_sub <- x[subset_indices[[k]]]

        message(sprintf("      subset %02d: fitting", k))

        fit_obj <- fit_single_chain(
          x = x_sub,
          init = init,
          sd_prop = sd_prop,
          temper = K,
          n_iter = n_iter,
          true_pars = true_pars,
          subset_index = k,
          K = K,
          store_data = store_data_with_fit
        )

        fit_obj$run_info$scenario_id <- scenario_id
        fit_obj$run_info$sim_id <- sim_id
        fit_obj$run_info$source_data_path <- file.path(sim_dir, "data.rds")
        fit_obj$run_info$source_meta_path <- file.path(sim_dir, "meta.rds")
        fit_obj$run_info$fit_type <- "subset"

        save_fit_bundle(
          fit_obj = fit_obj$draws,
          out_path = subset_out_path,
          run_info = fit_obj$run_info
        )
      }
    }
  }
}

message("")
message("Posterior fitting stage complete.")