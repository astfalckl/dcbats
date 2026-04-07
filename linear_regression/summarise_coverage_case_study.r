#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--root"), dest = "root", type = "character",
              default = "case_study_pilot",
              help = "Root directory containing scenario/replicate outputs [default %default]"),
  make_option(c("--n-reps"), dest = "n_reps", type = "integer",
              default = 3,
              help = "Number of replicates per scenario [default %default]"),
  make_option(c("--level"), dest = "level", type = "double",
              default = 0.95,
              help = "Credible interval level in (0,1) [default %default]"),
  make_option(c("--output-csv"), dest = "output_csv", type = "character",
              default = "case_study_pilot/coverage_summary_long.csv",
              help = "Output CSV for long-format coverage results [default %default]"),
  make_option(c("--output-table"), dest = "output_table", type = "character",
              default = "case_study_pilot/coverage_summary_table.csv",
              help = "Output CSV for paper-style summary table [default %default]")
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

coverage_indicators <- function(draws, beta_true, level = 0.95) {
  alpha <- 1 - level
  lower <- apply(draws, 2, quantile, probs = alpha / 2, names = FALSE, type = 7)
  upper <- apply(draws, 2, quantile, probs = 1 - alpha / 2, names = FALSE, type = 7)

  data.frame(
    coefficient = seq_along(beta_true),
    truth = as.numeric(beta_true),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    covered = as.integer((beta_true >= lower) & (beta_true <= upper))
  )
}

scenario_label <- function(x) {
  switch(
    x,
    iid = "i.i.d.",
    case1 = "Case I",
    case2 = "Case II",
    case3 = "Case III",
    x
  )
}

one_rep_summary <- function(root, scenario, rep, level) {
  rep_id <- sprintf("rep_%03d", rep)
  rep_dir <- file.path(root, scenario, rep_id)

  data_file <- file.path(rep_dir, "sim_data.rds")
  full_file <- file.path(rep_dir, "full_posterior_draws.rds")
  w5_file   <- file.path(rep_dir, "wasserstein_beta_K5.rds")
  w10_file  <- file.path(rep_dir, "wasserstein_beta_K10.rds")
  w20_file  <- file.path(rep_dir, "wasserstein_beta_K20.rds")

  check_file_exists(data_file, "Simulated data file")
  check_file_exists(full_file, "Full posterior file")
  check_file_exists(w5_file, "Wasserstein K=5 file")
  check_file_exists(w10_file, "Wasserstein K=10 file")
  check_file_exists(w20_file, "Wasserstein K=20 file")

  sim <- readRDS(data_file)
  full_object <- readRDS(full_file)
  w5_object <- readRDS(w5_file)
  w10_object <- readRDS(w10_file)
  w20_object <- readRDS(w20_file)

  if (!"beta_true" %in% names(sim)) {
    stop(sprintf("beta_true missing in %s", data_file), call. = FALSE)
  }

  beta_true <- as.numeric(sim$beta_true)
  beta_full <- extract_full_beta(full_object)
  beta_w5   <- extract_wasserstein_beta(w5_object)
  beta_w10  <- extract_wasserstein_beta(w10_object)
  beta_w20  <- extract_wasserstein_beta(w20_object)

  p <- length(beta_true)

  if (ncol(beta_full) != p || ncol(beta_w5) != p || ncol(beta_w10) != p || ncol(beta_w20) != p) {
    stop(sprintf("Dimension mismatch in %s", rep_dir), call. = FALSE)
  }

  full_cov <- coverage_indicators(beta_full, beta_true, level = level)
  w5_cov   <- coverage_indicators(beta_w5, beta_true, level = level)
  w10_cov  <- coverage_indicators(beta_w10, beta_true, level = level)
  w20_cov  <- coverage_indicators(beta_w20, beta_true, level = level)

  full_k5 <- transform(full_cov, scenario = scenario, replicate = rep, K = 5,  method = "Full")
  full_k10 <- transform(full_cov, scenario = scenario, replicate = rep, K = 10, method = "Full")
  full_k20 <- transform(full_cov, scenario = scenario, replicate = rep, K = 20, method = "Full")

  dc_k5  <- transform(w5_cov,  scenario = scenario, replicate = rep, K = 5,  method = "DC")
  dc_k10 <- transform(w10_cov, scenario = scenario, replicate = rep, K = 10, method = "DC")
  dc_k20 <- transform(w20_cov, scenario = scenario, replicate = rep, K = 20, method = "DC")

  rbind(dc_k5, full_k5, dc_k10, full_k10, dc_k20, full_k20)
}

paper_table_from_long <- function(long_df) {
  agg <- aggregate(
    covered ~ scenario + K + method,
    data = long_df,
    FUN = mean
  )

  agg$coverage_pct <- round(100 * agg$covered)
  agg$scenario_label <- vapply(agg$scenario, scenario_label, character(1))

  Ks <- c(5, 10, 20)
  scenarios <- c("iid", "case1", "case2", "case3")
  methods <- c("DC", "Full")

  out <- data.frame(
    K = paste0("K=", Ks),
    iid_DC = NA_integer_,
    iid_Full = NA_integer_,
    case1_DC = NA_integer_,
    case1_Full = NA_integer_,
    case2_DC = NA_integer_,
    case2_Full = NA_integer_,
    case3_DC = NA_integer_,
    case3_Full = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(Ks)) {
    Kval <- Ks[i]
    for (sc in scenarios) {
      for (m in methods) {
        hit <- agg$coverage_pct[agg$K == Kval & agg$scenario == sc & agg$method == m]
        if (length(hit) != 1L) {
          stop(sprintf("Missing or duplicated summary for scenario=%s, K=%d, method=%s", sc, Kval, m),
               call. = FALSE)
        }
        col_nm <- paste0(sc, "_", m)
        out[i, col_nm] <- hit
      }
    }
  }

  out
}

main <- function() {
  parser <- OptionParser(option_list = option_list)
  opt <- parse_args(parser)

  if (opt$n_reps < 1L) {
    stop("--n-reps must be positive.", call. = FALSE)
  }
  if (opt$level <= 0 || opt$level >= 1) {
    stop("--level must lie strictly between 0 and 1.", call. = FALSE)
  }

  scenarios <- c("iid", "case1", "case2", "case3")
  pieces <- vector("list", length(scenarios) * opt$n_reps)
  idx <- 1L

  for (sc in scenarios) {
    for (rep in seq_len(opt$n_reps)) {
      message(sprintf("Processing scenario=%s, rep=%03d", sc, rep))
      pieces[[idx]] <- one_rep_summary(
        root = opt$root,
        scenario = sc,
        rep = rep,
        level = opt$level
      )
      idx <- idx + 1L
    }
  }

  long_df <- do.call(rbind, pieces)
  long_df <- long_df[, c("scenario", "replicate", "K", "method",
                         "coefficient", "truth", "lower", "upper", "covered")]

  dir.create(dirname(opt$output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(long_df, file = opt$output_csv, row.names = FALSE)

  summary_table <- paper_table_from_long(long_df)
  dir.create(dirname(opt$output_table), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_table, file = opt$output_table, row.names = FALSE)

  cat("\nEmpirical coverage summary (%):\n")
  print(summary_table, row.names = FALSE)

  cat("\nPilot denominator per table entry:\n")
  cat(sprintf("n_reps x p = %d x %d = %d intervals per scenario/K/method\n",
              opt$n_reps,
              length(unique(long_df$coefficient)),
              opt$n_reps * length(unique(long_df$coefficient))))
}

main()