#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option(c("--data"), dest = "data", type = "character",
              default = "data/ar2_sim.rds",
              help = "Simulated data .rds from simulate_data.R [default %default]"),
  make_option(c("--full"), dest = "full", type = "character",
              default = "results/full_posterior_draws.rds",
              help = "Full posterior .rds from fit_posteriors.R [default %default]"),
  make_option(c("--w10"), dest = "w10", type = "character",
              default = "results/wasserstein_beta_K10.rds",
              help = "Wasserstein beta .rds for K = 10 [default %default]"),
  make_option(c("--w20"), dest = "w20", type = "character",
              default = "results/wasserstein_beta_K20.rds",
              help = "Wasserstein beta .rds for K = 20 [default %default]"),
  make_option(c("--output"), dest = "output", type = "character",
              default = "results/CI_ar2_errors_paper.pdf",
              help = "Output PDF path [default %default]"),
  make_option(c("--level"), dest = "level", type = "double",
              default = 0.95,
              help = "Credible interval level in (0,1) [default %default]"),
  make_option(c("--n-show"), dest = "n_show", type = "integer",
              default = 10,
              help = "Number of first coefficients to plot [default %default]"),
  make_option(c("--width"), dest = "width", type = "double",
              default = 10,
              help = "Figure width in inches [default %default]"),
  make_option(c("--height"), dest = "height", type = "double",
              default = 4.5,
              help = "Figure height in inches [default %default]"),
  make_option(c("--point-size"), dest = "point_size", type = "double",
              default = 1.8,
              help = "Point size for posterior mean markers [default %default]"),
  make_option(c("--truth-size"), dest = "truth_size", type = "double",
              default = 2.2,
              help = "Point size for true beta markers [default %default]")
)

check_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

extract_full_beta <- function(full_object) {
  if (!is.list(full_object) || !"draws" %in% names(full_object)) {
    stop("Full posterior object is malformed.", call. = FALSE)
  }
  if (!"beta" %in% names(full_object$draws)) {
    stop("Full posterior object does not contain draws$beta.", call. = FALSE)
  }
  as.matrix(full_object$draws$beta)
}

extract_wasserstein_beta <- function(w_object) {
  if (!is.list(w_object) || !"barycenter_draws" %in% names(w_object)) {
    stop("Wasserstein object is malformed.", call. = FALSE)
  }
  as.matrix(w_object$barycenter_draws)
}

# compute_interval_df <- function(draws, truth, method, level, idx) {
#   alpha <- 1 - level
#   lower_prob <- alpha / 2
#   upper_prob <- 1 - alpha / 2

#   lower <- apply(draws, 2, quantile, probs = lower_prob, names = FALSE, type = 7)
#   upper <- apply(draws, 2, quantile, probs = upper_prob, names = FALSE, type = 7)
#   center <- colMeans(draws)

#   data.frame(
#     index = idx,
#     truth = as.numeric(truth),
#     mean = as.numeric(center),
#     lower = as.numeric(lower),
#     upper = as.numeric(upper),
#     method = method,
#     stringsAsFactors = FALSE
#   )
# }

make_draws_long <- function(mat, method_label) {
  df <- as.data.frame(mat, check.names = FALSE)
  df$draw <- seq_len(nrow(df))

  out <- pivot_longer(
    df,
    cols = -draw,
    names_to = "index",
    values_to = "value"
  )

  out$index <- seq_len(ncol(mat))[match(out$index, colnames(mat))]
  out$method <- method_label
  out
}

# build_plot <- function(interval_df, truth_df, level, point_size, truth_size) {
#   dodge <- position_dodge(width = 0.6)

#   ggplot(interval_df, aes(x = factor(index), y = mean, ymin = lower, ymax = upper, shape = method)) +
#     geom_hline(yintercept = 0, linewidth = 0.25) +
#     geom_errorbar(position = dodge, width = 0.15, linewidth = 0.4) +
#     geom_point(position = dodge, size = point_size) +
#     geom_point(
#       data = truth_df,
#       aes(x = factor(index), y = truth),
#       inherit.aes = FALSE,
#       shape = 4,
#       size = truth_size,
#       stroke = 0.7
#     ) +
#     labs(
#       x = expression(beta~"index"),
#       y = expression(beta~"value"),
#       shape = NULL
#     ) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       plot.subtitle = element_text(hjust = 0.5),
#       legend.position = "right",
#       panel.grid.minor = element_blank()
#     )
# }

build_plot <- function(draws_df, truth_df, truth_size) {
  dodge <- position_dodge(width = 0.75)

  ggplot(draws_df, aes(x = factor(index), y = value, fill = method)) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_boxplot(
      position = dodge,
      width = 0.60,
      outlier.size = 0.35,
      linewidth = 0.35
    ) +
    geom_point(
      data = truth_df,
      aes(x = factor(index), y = truth),
      inherit.aes = FALSE,
      shape = 4,
      size = truth_size,
      stroke = 0.7
    ) +
    labs(
      x = expression(beta~"index"),
      y = expression(beta~"value"),
      fill = NULL
    ) +
    scale_fill_manual(values = c(
      "Wasserstein, K = 10" = "#4E79A7",  # blue
      "Wasserstein, K = 20" = "#F28E2B",  # orange
      "Full posterior"     = "#59A14F"    # green
    )) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

main <- function() {
  parser <- OptionParser(option_list = option_list)
  opt <- parse_args(parser)

  if (opt$level <= 0 || opt$level >= 1) {
    stop("--level must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (opt$n_show < 1) {
    stop("--n-show must be positive.", call. = FALSE)
  }
  if (opt$width <= 0 || opt$height <= 0) {
    stop("--width and --height must be positive.", call. = FALSE)
  }

  check_file_exists(opt$data, "Data file")
  check_file_exists(opt$full, "Full posterior file")
  check_file_exists(opt$w10, "Wasserstein K=10 file")
  check_file_exists(opt$w20, "Wasserstein K=20 file")

  sim <- readRDS(opt$data)
  full_object <- readRDS(opt$full)
  w10_object <- readRDS(opt$w10)
  w20_object <- readRDS(opt$w20)

  if (!"beta_true" %in% names(sim)) {
    stop("Simulated data object does not contain beta_true.", call. = FALSE)
  }

  beta_true <- as.numeric(sim$beta_true)
  beta_full <- extract_full_beta(full_object)
  beta_w10 <- extract_wasserstein_beta(w10_object)
  beta_w20 <- extract_wasserstein_beta(w20_object)

  p <- length(beta_true)
  if (ncol(beta_full) != p || ncol(beta_w10) != p || ncol(beta_w20) != p) {
    stop("Mismatch between beta_true length and posterior draw dimensions.", call. = FALSE)
  }

  keep_idx <- seq_len(min(opt$n_show, p))

  beta_true_sub <- beta_true[keep_idx]
  beta_full_sub <- beta_full[, keep_idx, drop = FALSE]
  beta_w10_sub  <- beta_w10[, keep_idx, drop = FALSE]
  beta_w20_sub  <- beta_w20[, keep_idx, drop = FALSE]

  # df_w10 <- compute_interval_df(
  #   draws = beta_w10_sub,
  #   truth = beta_true_sub,
  #   method = "Wasserstein, K = 10",
  #   level = opt$level,
  #   idx = keep_idx
  # )
  # df_w20 <- compute_interval_df(
  #   draws = beta_w20_sub,
  #   truth = beta_true_sub,
  #   method = "Wasserstein, K = 20",
  #   level = opt$level,
  #   idx = keep_idx
  # )
  # df_full <- compute_interval_df(
  #   draws = beta_full_sub,
  #   truth = beta_true_sub,
  #   method = "Full posterior",
  #   level = opt$level,
  #   idx = keep_idx
  # )

  draws_df <- bind_rows(
    make_draws_long(beta_w10_sub, "Wasserstein, K = 10"),
    make_draws_long(beta_w20_sub, "Wasserstein, K = 20"),
    make_draws_long(beta_full_sub, "Full posterior")
  )

  draws_df$method <- factor(
    draws_df$method,
    levels = c("Wasserstein, K = 10", "Wasserstein, K = 20", "Full posterior")
  )

  # interval_df <- rbind(df_w10, df_w20, df_full)
  # interval_df$method <- factor(
  #   interval_df$method,
  #   levels = c("Wasserstein, K = 10", "Wasserstein, K = 20", "Full posterior")
  # )

  truth_df <- data.frame(
    index = keep_idx,
    truth = beta_true_sub
  )

  # fig <- build_plot(
  #   interval_df = interval_df,
  #   truth_df = truth_df,
  #   level = opt$level,
  #   point_size = opt$point_size,
  #   truth_size = opt$truth_size
  # )

  fig <- build_plot(
    draws_df = draws_df,
    truth_df = truth_df,
    truth_size = opt$truth_size
  )

  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = opt$output, plot = fig, width = opt$width, height = opt$height)

  message(sprintf("Saved figure to %s", opt$output))
}

main()