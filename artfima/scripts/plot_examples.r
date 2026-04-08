# ------------------------------------------------------------------------------
# Plot Wasserstein barycenter distributions for one simulation
# ------------------------------------------------------------------------------

rm(list = ls())

source("R/helpers_artfima.R")
load_required_packages()
library(ggplot2)

source("R/wasserstein_utils.R")

fit_root <- file.path("results", "fits")

parameters <- c("phi")
probs <- seq(0.001, 0.999, by = 0.001)
K_values <- c(5L, 10L, 20L)

sim_id_target <- "012"
scenario_ids <- c("lambda0005_d01_phi01",
                  "lambda0005_d01_phi09",
                  "lambda0005_d01_phi099")

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

# In 1D, the Wasserstein barycenter is the pointwise mean of quantile functions.
compute_barycenter_quantiles <- function(draws_list, param, probs) {
  lapply(draws_list, function(d) {
    compute_quantile_df(d, param, probs)
  }) |>
    bind_rows(.id = "replicate_id") |>
    group_by(p) |>
    summarise(value = mean(value), .groups = "drop")
}

# Convert a quantile curve into pseudo-draws for density plotting.
# This samples U ~ Unif(0,1) and returns Q(U) by interpolation.
quantile_df_to_draws <- function(qdf, n = 5000) {
  u <- runif(n)
  approx(
    x = qdf$p,
    y = qdf$value,
    xout = u,
    ties = mean,
    rule = 2
  )$y
}

all_plot_data <- list()

for (scenario_id in scenario_ids) {

  sim_dir <- file.path(fit_root, scenario_id, sim_id_target)

  full_path <- file.path(sim_dir, "full.rds")
  full_draws <- read_draws(full_path)

  for (param in parameters) {

    # Full posterior: use actual draws
    all_plot_data[[length(all_plot_data) + 1L]] <- tibble(
      scenario_id = scenario_id,
      sim_id = sim_id_target,
      parameter = param,
      source = "Full Post.",
      value = full_draws[, param]
    )

    # Wasserstein barycenters for each K
    for (K in K_values) {

      subset_dir <- file.path(sim_dir, sprintf("K%02d", K))
      files <- list.files(subset_dir, full.names = TRUE)

      draws_list <- lapply(files, read_draws)

      bary_q <- compute_barycenter_quantiles(draws_list, param, probs)
      bary_draws <- quantile_df_to_draws(bary_q, n = 5000)

      all_plot_data[[length(all_plot_data) + 1L]] <- tibble(
        scenario_id = scenario_id,
        sim_id = sim_id_target,
        parameter = param,
        source = sprintf("K=%d", K),
        value = bary_draws
      )
    }
  }
}

plot_data <- bind_rows(all_plot_data) |>
  mutate(
    phi_label = case_when(
      grepl("phi01$", scenario_id)  ~ "phi == 0.1",
      grepl("phi09$", scenario_id)  ~ "phi == 0.9",
      grepl("phi099$", scenario_id) ~ "phi == 0.99",
      TRUE ~ scenario_id
    ),
    source = factor(
      source,
      levels = c("Full Post.", "K=5", "K=10", "K=20")
    )
  )

# Cleaner facet labels
facet_labs <- c(
  "lambda0005_d01_phi01"  = "phi == 0.1",
  "lambda0005_d01_phi09"  = "phi == 0.9",
  "lambda0005_d01_phi099" = "phi == 0.99"
)

true_phi <- tibble::tibble(
  scenario_id = c(
    "lambda0005_d01_phi01",
    "lambda0005_d01_phi09",
    "lambda0005_d01_phi099"
  ),
  phi_true = c(0.1, 0.9, 0.99)
)

cols <- c(
  "Full Post."           = "#1F2933",  # charcoal
  "K=5" = "#4C78A8",  # muted blue
  "K=10" = "#F58518",  # muted orange
  "K=20" = "#54A24B"   # muted green
)

p <- ggplot(plot_data, aes(x = value, colour = source, fill = source)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.5) +
  facet_wrap(
    ~ scenario_id,
    scales = "free",
    labeller = as_labeller(facet_labs, label_parsed)
  ) +
  geom_vline(
    data = true_phi,
    aes(xintercept = phi_true),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = expression(phi),
    y = "Density",
    colour = NULL,
    fill = NULL,
  ) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.line = element_line(colour = "grey50"),
    axis.ticks = element_line(colour = "grey50")
  )

print(p)

ggsave(
  filename = file.path(
    "results", "summaries",
    paste0("phi_barycenter_kde_sim", sim_id_target, ".pdf")
  ),
  plot = p,
  width = 8,
  height = 2.5,
  dpi = 300
)
