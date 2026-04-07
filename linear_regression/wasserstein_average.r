#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

# Compute coefficientwise equal-weight Wasserstein barycenter draws for beta
# from subset posterior draws.
#
# For each coefficient j and rank i, the barycenter draw is
#   beta_star_(i,j) = (1/K) sum_{k=1}^K beta_(i,j)^(k)
# where beta_(i,j)^(k) is the i-th order statistic of the subset-k draws for
# coefficient j.
#
# This script expects an .rds file produced by fit_posteriors.R, specifically
# one of the subset_posteriors_K*.rds files.

option_list <- list(
  make_option(c("--input"), dest = "input", type = "character",
              default = "results/subset_posteriors_K10.rds",
              help = "Subset posterior .rds from fit_posteriors.R [default %default]"),
  make_option(c("--output"), dest = "output", type = "character",
              default = "results/wasserstein_beta_K10.rds",
              help = "Output .rds file for barycenter draws [default %default]"),
  make_option(c("--parameter"), dest = "parameter", type = "character",
              default = "beta",
              help = "Parameter to aggregate; currently only 'beta' is supported [default %default]"),
  make_option(c("--probs"), dest = "probs", type = "character",
              default = "0.025,0.5,0.975",
              help = "Comma-separated posterior summary probabilities [default %default]"),
  make_option(c("--check-dimnames"), dest = "check_dimnames", action = "store_true",
              default = TRUE,
              help = "Check column names across subset draw matrices [default %default]"),
  make_option(c("--no-check-dimnames"), action = "store_false", dest = "check_dimnames",
              help = "Do not check column names across subset draw matrices")
)

parse_probs <- function(prob_string) {
  probs <- trimws(strsplit(prob_string, split = ",", fixed = TRUE)[[1]])
  probs <- as.numeric(probs[nzchar(probs)])
  if (length(probs) < 1L || any(is.na(probs)) || any(probs <= 0) || any(probs >= 1)) {
    stop("--probs must be a comma-separated list of numbers strictly between 0 and 1.", call. = FALSE)
  }
  unique(probs)
}

validate_subset_object <- function(x) {
  needed <- c("K", "blocks", "fits", "meta")
  if (!all(needed %in% names(x))) {
    stop("Input object is not a valid subset posterior object.", call. = FALSE)
  }
  if (!is.list(x$fits) || length(x$fits) != x$K) {
    stop("Input object has invalid subset fit structure.", call. = FALSE)
  }
}

extract_parameter_draws <- function(subset_object, parameter = "beta", check_dimnames = TRUE) {
  if (parameter != "beta") {
    stop("Only parameter = 'beta' is currently supported.", call. = FALSE)
  }

  draws_list <- lapply(subset_object$fits, function(f) {
    if (!is.list(f) || !"draws" %in% names(f) || !parameter %in% names(f$draws)) {
      stop("A subset fit does not contain the requested parameter draws.", call. = FALSE)
    }
    as.matrix(f$draws[[parameter]])
  })

  n_rows <- vapply(draws_list, nrow, integer(1))
  n_cols <- vapply(draws_list, ncol, integer(1))

  if (length(unique(n_rows)) != 1L) {
    stop("Subset draw matrices have differing numbers of rows.", call. = FALSE)
  }
  if (length(unique(n_cols)) != 1L) {
    stop("Subset draw matrices have differing numbers of columns.", call. = FALSE)
  }

  if (isTRUE(check_dimnames)) {
    ref_names <- colnames(draws_list[[1]])
    for (k in seq_along(draws_list)) {
      this_names <- colnames(draws_list[[k]])
      if (!identical(ref_names, this_names)) {
        stop(sprintf("Column names differ in subset %d.", k), call. = FALSE)
      }
    }
  }

  draws_list
}

wasserstein_average_matrix <- function(draws_list) {
  K <- length(draws_list)
  n_draws <- nrow(draws_list[[1]])
  p <- ncol(draws_list[[1]])

  out <- matrix(NA_real_, nrow = n_draws, ncol = p)

  for (j in seq_len(p)) {
    sorted_by_subset <- vapply(
      draws_list,
      function(mat) sort(mat[, j]),
      numeric(n_draws)
    )
    out[, j] <- rowMeans(sorted_by_subset)
  }

  colnames(out) <- colnames(draws_list[[1]])
  out
}

summarise_draws <- function(draws, probs) {
  qmat <- t(apply(draws, 2, quantile, probs = probs, names = FALSE, type = 7))
  out <- data.frame(
    parameter = colnames(draws),
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(probs)) {
    nm <- paste0("q", formatC(probs[i] * 100, format = "f", digits = 1))
    out[[nm]] <- qmat[, i]
  }

  out
}

save_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = path)
}

main <- function() {
  parser <- OptionParser(option_list = option_list)
  opt <- parse_args(parser)

  if (!file.exists(opt$input)) {
    stop(sprintf("Input file not found: %s", opt$input), call. = FALSE)
  }
  if (!identical(opt$parameter, "beta")) {
    stop("Currently only --parameter beta is supported.", call. = FALSE)
  }

  probs <- parse_probs(opt$probs)
  subset_object <- readRDS(opt$input)
  validate_subset_object(subset_object)

  draws_list <- extract_parameter_draws(
    subset_object = subset_object,
    parameter = opt$parameter,
    check_dimnames = opt$check_dimnames
  )

  beta_barycenter <- wasserstein_average_matrix(draws_list)
  summary_df <- summarise_draws(beta_barycenter, probs = probs)

  out <- list(
    parameter = opt$parameter,
    K = subset_object$K,
    barycenter_draws = beta_barycenter,
    summary = summary_df,
    meta = list(
      source_file = normalizePath(opt$input, winslash = "/", mustWork = FALSE),
      probs = probs,
      n_subsets = length(draws_list),
      n_draws = nrow(beta_barycenter),
      p = ncol(beta_barycenter),
      method = "coefficientwise equal-weight Wasserstein barycenter via order-statistic averaging"
    )
  )

  save_object(out, opt$output)

  message(sprintf("Saved Wasserstein-averaged draws to %s", opt$output))
  message(sprintf("K = %d, n_draws = %d, p = %d", out$K, nrow(beta_barycenter), ncol(beta_barycenter)))
}

main()
