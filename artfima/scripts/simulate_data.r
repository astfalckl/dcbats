# ------------------------------------------------------------------------------
# 01_simulate_data.R
#
# Simulate ARTFIMA datasets for the reproducibility pipeline.
# This script:
#   1. loads shared helpers
#   2. defines the simulation design
#   3. simulates one full dataset per scenario/replicate
#   4. saves data + metadata only
#
# It does NOT fit posteriors.
# ------------------------------------------------------------------------------

rm(list = ls())

# ------------------------------------------------------------------------------
# Load shared helpers
# ------------------------------------------------------------------------------

source("R/helpers_artfima.R")
load_required_packages()
source_artfima_code()

source("R/fast_acf_from_psd.R")

# ------------------------------------------------------------------------------
# User-configurable settings
# ------------------------------------------------------------------------------

seed_global <- 20260407

n_data <- 10000
n_simulations <- 100

# Root for simulated data only
simulation_root <- file.path("data", "simulated")

# Whether to overwrite an existing data.rds
overwrite_existing <- FALSE

# Use default grid from helper unless you want to replace it here
param_grid <- default_param_grid()

# ------------------------------------------------------------------------------
# Directory setup
# ------------------------------------------------------------------------------

ensure_project_subdirs(c(
  "data",
  simulation_root
))

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

simulate_artfima_dataset <- function(n, d, lambda, phi, sigma2 = 1) {
  acf_design <- bigdog_artfima_acf(
    n = n - 1L,
    d = d,
    lambda = lambda,
    phi = phi
  )

  x <- rnormtz(
    n = 1,
    acf_design,
    fft = FALSE
  )

  list(
    x = x,
    acf = acf_design,
    sigma2 = sigma2
  )
}

write_simulation_bundle <- function(
  out_dir,
  data_vec,
  scenario,
  sim_index,
  seed,
  n_data,
  acf_vec = NULL,
  overwrite = FALSE
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  data_path <- file.path(out_dir, "data.rds")
  meta_path <- file.path(out_dir, "meta.rds")

  if (!overwrite && file.exists(data_path) && file.exists(meta_path)) {
    message(sprintf("Skipping existing simulation: %s", out_dir))
    return(invisible(FALSE))
  }

  metadata <- list(
    scenario = scenario,
    sim_index = sim_index,
    seed = seed,
    n_data = n_data,
    generated_at = Sys.time(),
    generator = "bigdog_artfima_acf + rnormtz",
    files = list(
      data = data_path,
      meta = meta_path
    )
  )

  saveRDS(data_vec, data_path)
  saveRDS(metadata, meta_path)

  # Optional: save the design acf for auditing/debugging
  if (!is.null(acf_vec)) {
    saveRDS(acf_vec, file.path(out_dir, "acf_design.rds"))
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

set.seed(seed_global)

n_scenarios <- nrow(param_grid)

for (i in seq_len(n_scenarios)) {
  scenario <- param_grid[i, , drop = FALSE]

  scenario_id <- make_scenario_id(
    lambda = scenario$lambda,
    d = scenario$d,
    phi = scenario$phi
  )

  scenario_dir <- file.path(simulation_root, scenario_id)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)

  message("")
  message(sprintf(
    "Scenario %d/%d: %s (lambda = %s, d = %s, phi = %s)",
    i,
    n_scenarios,
    scenario_id,
    scenario$lambda,
    scenario$d,
    scenario$phi
  ))

  for (j in seq_len(n_simulations)) {
    sim_id <- make_sim_id(j)
    sim_dir <- file.path(scenario_dir, sim_id)

    sim_seed <- seed_global + 100000L * i + j
    set.seed(sim_seed)

    message(sprintf("  Sim %s", sim_id))

    sim_obj <- simulate_artfima_dataset(
      n = n_data,
      d = scenario$d,
      lambda = scenario$lambda,
      phi = scenario$phi,
      sigma2 = 1
    )

    write_simulation_bundle(
      out_dir = sim_dir,
      data_vec = sim_obj$x,
      scenario = as.list(scenario),
      sim_index = j,
      seed = sim_seed,
      n_data = n_data,
      acf_vec = sim_obj$acf,
      overwrite = overwrite_existing
    )
  }
}

message("")
message("Simulation stage complete.")
