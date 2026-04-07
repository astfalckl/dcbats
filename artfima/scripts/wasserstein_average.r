# ------------------------------------------------------------------------------
# 03_wasserstein_average.R
#
# Compute posterior quantiles and Wasserstein barycenter approximations.
#
# Input: results/fits/
# Output: results/summaries/posterior_quantiles_<param>.rds
# ------------------------------------------------------------------------------

rm(list = ls())

# ------------------------------------------------------------------------------
# Load helpers
# ------------------------------------------------------------------------------

source("R/helpers_artfima.R")
load_required_packages()

source("R/wasserstein_utils.R")

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

fit_root <- file.path("results", "fits")
summary_root <- file.path("results", "summaries")

ensure_project_subdirs(summary_root)

parameters <- c("phi")

probs <- seq(0.01, 0.99, by = 0.01)

K_values <- c(5L, 10L, 20L)

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

read_draws <- function(path) {
  obj <- readRDS(path)
  obj$draws
}

compute_quantile_df <- function(draws, param, probs) {
  tibble(
    p = probs,
    value = compute_percentiles(draws[, param], probs = probs)
  )
}

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

scenario_dirs <- list.dirs(fit_root, recursive = FALSE)

# scenario_dirs <- scenario_dirs[1:6]

all_results <- list()

for (scenario_dir in scenario_dirs) {

  scenario_id <- basename(scenario_dir)
  sim_dirs <- list.dirs(scenario_dir, recursive = FALSE)

  message(sprintf("Scenario: %s", scenario_id))

  for (sim_dir in sim_dirs) {

    sim_id <- basename(sim_dir)

    message(sprintf("  Simulation: %s", sim_id))

    # --------------------------
    # FULL posterior
    # --------------------------

    full_path <- file.path(sim_dir, "full.rds")
    full_draws <- read_draws(full_path)

    # --------------------------
    # SUBSET posteriors
    # --------------------------

    subset_draws <- list()

    for (K in K_values) {

      subset_dir <- file.path(sim_dir, sprintf("K%02d", K))
      files <- list.files(subset_dir, full.names = TRUE)

      subset_draws[[paste0("K", K)]] <- lapply(files, read_draws)
    }

    # --------------------------
    # Parameter loop
    # --------------------------

    for (param in parameters) {

      full_q <- compute_quantile_df(full_draws, param, probs) %>%
        rename(full = value)

      subset_q <- lapply(names(subset_draws), function(kname) {

        draws_list <- subset_draws[[kname]]

        df <- lapply(draws_list, function(d) {
          compute_quantile_df(d, param, probs)
        }) %>%
          bind_rows() %>%
          group_by(p) %>%
          summarise(value = mean(value), .groups = "drop")

        df %>% rename(!!kname := value)
      }) %>%
        reduce(left_join, by = "p")

      result <- full_q %>%
        left_join(subset_q, by = "p") %>%
        mutate(
          scenario_id = scenario_id,
          sim_id = sim_id,
          parameter = param
        ) %>%
        select(scenario_id, sim_id, parameter, p, everything())

      all_results[[length(all_results) + 1L]] <- result
    }
  }
}

posterior_quantiles <- bind_rows(all_results)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

out_path <- file.path(summary_root, "posterior_quantiles.rds")
saveRDS(posterior_quantiles, out_path)

message("Saved: ", out_path)
