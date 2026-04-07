# ------------------------------------------------------------------------------
# 04_aggregate_results.R
#
# Aggregate Wasserstein discrepancies across simulation replicates.
#
# Input:
#   results/summaries/posterior_quantiles.rds
#
# Output:
#   results/summaries/wasserstein_by_sim.rds
#   results/summaries/wasserstein_summary.rds
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

summary_root <- file.path("results", "summaries")
quantile_path <- file.path(summary_root, "posterior_quantiles.rds")

if (!file.exists(quantile_path)) {
  stop(sprintf("Missing file: %s", quantile_path), call. = FALSE)
}

distance_method <- "riemann"

# ------------------------------------------------------------------------------
# Read quantile summaries
# ------------------------------------------------------------------------------

posterior_quantiles <- readRDS(quantile_path)

required_cols <- c("scenario_id", "sim_id", "parameter", "p", "full", "K5", "K10", "K20")
missing_cols <- setdiff(required_cols, names(posterior_quantiles))

if (length(missing_cols) > 0L) {
  stop(
    sprintf(
      "posterior_quantiles is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Compute Wasserstein discrepancy for each simulation replicate
# ------------------------------------------------------------------------------

wasserstein_by_sim <- posterior_quantiles %>%
  dplyr::group_by(scenario_id, sim_id, parameter) %>%
  dplyr::summarise(
    wass_K5 = W1_rel(full, K5, method = distance_method),
    wass_K10 = W1_rel(full, K10, method = distance_method),
    wass_K20 = W1_rel(full, K20, method = distance_method),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Aggregate across replicates
# ------------------------------------------------------------------------------

wasserstein_summary <- wasserstein_by_sim %>%
  dplyr::group_by(scenario_id, parameter) %>%
  dplyr::summarise(
    n_sim = dplyr::n(),

    mean_wass_K5 = mean(wass_K5),
    mean_wass_K10 = mean(wass_K10),
    mean_wass_K20 = mean(wass_K20),

    # median_wass_K5 = median(wass_K5),
    # median_wass_K10 = median(wass_K10),
    # median_wass_K20 = median(wass_K20),

    # sd_wass_K5 = stats::sd(wass_K5),
    # sd_wass_K10 = stats::sd(wass_K10),
    # sd_wass_K20 = stats::sd(wass_K20),

    # se_mean_wass_K5 = sd_wass_K5 / sqrt(n_sim),
    # se_mean_wass_K10 = sd_wass_K10 / sqrt(n_sim),
    # se_mean_wass_K20 = sd_wass_K20 / sqrt(n_sim),

    # min_wass_K5 = min(wass_K5),
    # min_wass_K10 = min(wass_K10),
    # min_wass_K20 = min(wass_K20),

    # max_wass_K5 = max(wass_K5),
    # max_wass_K10 = max(wass_K10),
    # max_wass_K20 = max(wass_K20),

    .groups = "drop"
  ) %>%
  dplyr::arrange(parameter, scenario_id)

# ------------------------------------------------------------------------------
# Optional: reshape to long form for plotting/tables
# ------------------------------------------------------------------------------

wasserstein_summary_long <- wasserstein_by_sim %>%
  tidyr::pivot_longer(
    cols = c(wass_K5, wass_K10, wass_K20),
    names_to = "K",
    values_to = "wasserstein_rel"
  ) %>%
  dplyr::mutate(
    K = dplyr::recode(K,
      wass_K5 = "K5",
      wass_K10 = "K10",
      wass_K20 = "K20"
    )
  )

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

saveRDS(
  wasserstein_by_sim,
  file.path(summary_root, "wasserstein_by_sim.rds")
)

saveRDS(
  wasserstein_summary,
  file.path(summary_root, "wasserstein_summary.rds")
)

saveRDS(
  wasserstein_summary_long,
  file.path(summary_root, "wasserstein_summary_long.rds")
)

message("Saved:")
message("  ", file.path(summary_root, "wasserstein_by_sim.rds"))
message("  ", file.path(summary_root, "wasserstein_summary.rds"))
message("  ", file.path(summary_root, "wasserstein_summary_long.rds"))