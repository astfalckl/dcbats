#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

# Simulate data for the linear regression model with AR(2) errors:
#   y_t = alpha + z_t^T beta + eps_t
#   eps_t = phi1 * eps_{t-1} + phi2 * eps_{t-2} + xi_t
#   xi_t ~ N(0, sigma2)
#
# The manuscript states independent N(0, sigma^2) initial conditions for
# eps_0 and eps_{-1}; in 1-based indexing below these correspond to eps[2]
# and eps[1], respectively.
#
# This script writes a single .rds file containing the simulated data and the
# true parameter values.

option_list <- list(
  make_option("--T", type = "integer", default = 100000,
              help = "Number of observations [default %default]"),
  make_option("--p", type = "integer", default = 50,
              help = "Number of covariates [default %default]"),
  make_option("--alpha", type = "double", default = 0,
              help = "Intercept alpha [default %default]"),
  make_option("--phi1", type = "double", default = 0.4,
              help = "AR(2) coefficient phi1 [default %default]"),
  make_option("--phi2", type = "double", default = -0.6,
              help = "AR(2) coefficient phi2 [default %default]"),
  make_option("--sigma2", type = "double", default = 1,
              help = "Innovation variance sigma^2 [default %default]"),
  make_option(c("--beta-sd"), dest = "beta_sd", type = "double", default = 0.1,
              help = "Standard deviation used to simulate beta_j [default %default]"),
  make_option(c("--z-dist"), dest = "z_dist", type = "character", default = "normal",
              help = "Covariate distribution: normal or uniform [default %default]"),
  make_option("--seed", type = "integer", default = 20260402,
              help = "Random seed [default %default]"),
  make_option("--output", type = "character", default = "data/sim_data.rds",
              help = "Output .rds path [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

validate_inputs <- function(opt) {
  if (opt$T < 3L) stop("T must be at least 3.", call. = FALSE)
  if (opt$p < 1L) stop("p must be positive.", call. = FALSE)
  if (opt$sigma2 <= 0) stop("sigma2 must be positive.", call. = FALSE)
  if (opt$beta_sd <= 0) stop("beta-sd must be positive.", call. = FALSE)
  if (!opt$z_dist %in% c("normal", "uniform")) {
    stop("z-dist must be one of: normal, uniform.", call. = FALSE)
  }

  # AR(2) stationarity conditions for 1 - phi1 z - phi2 z^2 = 0
  phi1 <- opt$phi1
  phi2 <- opt$phi2
  stationary <- (abs(phi2) < 1) && (phi1 + phi2 < 1) && (phi2 - phi1 < 1)
  if (!stationary) {
    warning(
      sprintf(
        paste0(
          "The supplied AR(2) coefficients may be nonstationary: ",
          "phi1 = %.4f, phi2 = %.4f. Proceeding anyway."
        ),
        phi1, phi2
      ),
      call. = FALSE
    )
  }
}

simulate_covariates <- function(T, p, dist) {
  if (dist == "normal") {
    matrix(rnorm(T * p), nrow = T, ncol = p)
  } else if (dist == "uniform") {
    matrix(runif(T * p, min = -1, max = 1), nrow = T, ncol = p)
  } else {
    stop("Unsupported covariate distribution.", call. = FALSE)
  }
}

simulate_ar2_errors <- function(T, phi1, phi2, sigma2) {
  eps <- numeric(T + 2L)

  # eps[1] = eps_{-1}, eps[2] = eps_0
  eps[1] <- rnorm(1L, mean = 0, sd = sqrt(sigma2))
  eps[2] <- rnorm(1L, mean = 0, sd = sqrt(sigma2))

  for (t in seq_len(T)) {
    idx <- t + 2L
    xi_t <- rnorm(1L, mean = 0, sd = sqrt(sigma2))
    eps[idx] <- phi1 * eps[idx - 1L] + phi2 * eps[idx - 2L] + xi_t
  }

  # return eps_1, ..., eps_T
  eps[3:(T + 2L)]
}

main <- function(opt) {
  validate_inputs(opt)
  set.seed(opt$seed)

  T <- opt$T
  p <- opt$p
  alpha <- opt$alpha
  phi1 <- opt$phi1
  phi2 <- opt$phi2
  sigma2 <- opt$sigma2
  beta_sd <- opt$beta_sd

  Z <- simulate_covariates(T = T, p = p, dist = opt$z_dist)
  beta <- rnorm(p, mean = 0, sd = beta_sd)
  eps <- simulate_ar2_errors(T = T, phi1 = phi1, phi2 = phi2, sigma2 = sigma2)
  y <- as.numeric(alpha + Z %*% beta + eps)

  out <- list(
    y = y,
    Z = Z,
    beta_true = beta,
    alpha_true = alpha,
    phi_true = c(phi1 = phi1, phi2 = phi2),
    sigma2_true = sigma2,
    eps = eps,
    meta = list(
      seed = opt$seed,
      T = T,
      p = p,
      z_dist = opt$z_dist,
      beta_sd = beta_sd,
      model = "linear regression with AR(2) errors"
    )
  )

  out_dir <- dirname(opt$output)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  saveRDS(out, file = opt$output)

  cat(sprintf("Wrote simulated data to %s\n", opt$output))
  cat(sprintf("T = %d, p = %d, seed = %d\n", T, p, opt$seed))
  cat(sprintf("alpha = %.4f, phi = (%.4f, %.4f), sigma2 = %.4f\n",
              alpha, phi1, phi2, sigma2))
}

main(opt)
