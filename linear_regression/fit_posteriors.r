#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(cmdstanr)
  library(posterior)
})

# Fit subset posteriors and the full posterior for the linear regression model
# with AR(2) errors. The same Stan model is used for both:
#   - full posterior: temper = 1
#   - subset posterior for K blocks: temper = K
#
# Input:
#   an .rds file produced by simulate_data.R
#
# Output:
#   - full posterior draws for beta and selected scalar parameters
#   - subset posterior draws for beta and selected scalar parameters
#   - metadata describing the split and MCMC settings
#
# Expected Stan data variables:
#   T, p, Z, y, temper

option_list <- list(
  make_option(c("--data"), dest = "data", type = "character", default = "data/sim_data.rds",
              help = "Input .rds file from simulate_data.R [default %default]"),
  make_option(c("--stan-file"), dest = "stan_file", type = "character", default = "model_ar2_errors.stan",
              help = "Stan model file [default %default]"),
  make_option(c("--output-dir"), dest = "output_dir", type = "character", default = "results",
              help = "Directory for posterior draw files [default %default]"),
  make_option(c("--K"), dest = "K", type = "character", default = "10,20",
              help = "Comma-separated subset counts, e.g. 10,20 [default %default]"),
  make_option(c("--chains"), dest = "chains", type = "integer", default = 4,
              help = "Number of MCMC chains [default %default]"),
  make_option(c("--parallel-chains"), dest = "parallel_chains", type = "integer", default = 4,
              help = "Number of parallel chains [default %default]"),
  make_option(c("--iter-warmup"), dest = "iter_warmup", type = "integer", default = 1000,
              help = "Warmup iterations per chain [default %default]"),
  make_option(c("--iter-sampling"), dest = "iter_sampling", type = "integer", default = 1000,
              help = "Post-warmup iterations per chain [default %default]"),
  make_option(c("--seed"), dest = "seed", type = "integer", default = 20260407,
              help = "Base random seed [default %default]"),
  make_option(c("--adapt-delta"), dest = "adapt_delta", type = "double", default = 0.9,
              help = "NUTS adapt_delta [default %default]"),
  make_option(c("--max-treedepth"), dest = "max_treedepth", type = "integer", default = 12,
              help = "NUTS max_treedepth [default %default]"),
  make_option(c("--refresh"), dest = "refresh", type = "integer", default = 100,
              help = "Stan sampling progress frequency [default %default]"),
  make_option(c("--save-full"), dest = "save_full", action = "store_true", default = TRUE,
              help = "Save full posterior draws [default %default]"),
  make_option("--no-save-full", action = "store_false", dest = "save_full",
              help = "Do not save full posterior draws"),
  make_option(c("--save-subsets"), dest = "save_subsets", action = "store_true", default = TRUE,
              help = "Save subset posterior draws [default %default]"),
  make_option("--no-save-subsets", action = "store_false", dest = "save_subsets",
              help = "Do not save subset posterior draws")
)

parse_K_values <- function(K_string) {
  K_vals <- trimws(strsplit(K_string, split = ",", fixed = TRUE)[[1]])
  K_vals <- as.integer(K_vals[nzchar(K_vals)])
  if (length(K_vals) < 1L || any(is.na(K_vals)) || any(K_vals < 1L)) {
    stop("--K must be a comma-separated list of positive integers.", call. = FALSE)
  }
  unique(K_vals)
}

# make_blocks <- function(T, K) {
#   # Nearly equal contiguous blocks that partition 1:T.
#   starts <- floor(seq.int(from = 0, to = T - 1, length.out = K)) + 1L
#   ends <- c(starts[-1L] - 1L, T)
#   sizes <- ends - starts + 1L

#   out <- data.frame(
#     block = seq_len(K),
#     start = starts,
#     end = ends,
#     size = sizes
#   )

#   if (sum(out$size) != T) {
#     stop("Block construction failed: sizes do not sum to T.", call. = FALSE)
#   }

#   out
# }

make_blocks <- function(T, K) {
  cuts <- floor(seq(0, T, length.out = K + 1L))
  starts <- cuts[-length(cuts)] + 1L
  ends <- cuts[-1L]
  sizes <- ends - starts + 1L

  data.frame(
    block = seq_len(K),
    start = starts,
    end = ends,
    size = sizes
  )
}

make_stan_data <- function(y, Z, temper) {
  list(
    T = length(y),
    p = ncol(Z),
    Z = unname(Z),
    y = as.numeric(y),
    temper = as.numeric(temper)
  )
}

extract_draws <- function(fit) {
  draws_df <- posterior::as_draws_df(fit$draws())

  beta_cols <- grep("^beta\\[", names(draws_df), value = TRUE)
  beta_mat <- as.matrix(draws_df[, beta_cols, drop = FALSE])

  scalar_candidates <- c("alpha", "phi1", "phi2", "sigma2")
  scalar_cols <- intersect(scalar_candidates, names(draws_df))
  scalar_df <- as.data.frame(draws_df[, scalar_cols, drop = FALSE])

  list(
    beta = beta_mat,
    scalars = scalar_df,
    diagnostics = fit$diagnostic_summary()
  )
}

fit_one_dataset <- function(model, y, Z, temper, seed, chains, parallel_chains,
                            iter_warmup, iter_sampling, adapt_delta,
                            max_treedepth, refresh) {
  stan_data <- make_stan_data(y = y, Z = Z, temper = temper)

  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    show_messages = TRUE
  )

  extract_draws(fit)
}

save_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = path)
}

main <- function() {
  parser <- OptionParser(option_list = option_list)
  opt <- parse_args(parser)

  K_values <- parse_K_values(opt$K)

  if (!file.exists(opt$data)) {
    stop(sprintf("Data file not found: %s", opt$data), call. = FALSE)
  }
  if (!file.exists(opt$stan_file)) {
    stop(sprintf("Stan file not found: %s", opt$stan_file), call. = FALSE)
  }
  if (opt$chains < 1L || opt$parallel_chains < 1L) {
    stop("chains and parallel-chains must be positive integers.", call. = FALSE)
  }
  if (opt$iter_warmup < 0L || opt$iter_sampling < 1L) {
    stop("iter-warmup must be nonnegative and iter-sampling must be positive.", call. = FALSE)
  }

  sim <- readRDS(opt$data)

  required_names <- c("y", "Z")
  if (!all(required_names %in% names(sim))) {
    stop("Input .rds must contain at least 'y' and 'Z'.", call. = FALSE)
  }

  y <- sim$y
  Z <- sim$Z
  T_total <- length(y)
  p <- ncol(Z)

  if (!is.numeric(y) || !is.matrix(Z)) {
    stop("Input data must have numeric vector y and numeric matrix Z.", call. = FALSE)
  }
  if (nrow(Z) != T_total) {
    stop("nrow(Z) must equal length(y).", call. = FALSE)
  }

  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Compiling Stan model...")
  model <- cmdstan_model(opt$stan_file)

  message("Fitting full posterior (temper = 1)...")
  full_fit <- fit_one_dataset(
    model = model,
    y = y,
    Z = Z,
    temper = 1,
    seed = opt$seed,
    chains = opt$chains,
    parallel_chains = opt$parallel_chains,
    iter_warmup = opt$iter_warmup,
    iter_sampling = opt$iter_sampling,
    adapt_delta = opt$adapt_delta,
    max_treedepth = opt$max_treedepth,
    refresh = opt$refresh
  )

  full_object <- list(
    draws = full_fit,
    meta = list(
      data_file = normalizePath(opt$data, winslash = "/", mustWork = FALSE),
      stan_file = normalizePath(opt$stan_file, winslash = "/", mustWork = FALSE),
      T = T_total,
      p = p,
      temper = 1,
      chains = opt$chains,
      parallel_chains = opt$parallel_chains,
      iter_warmup = opt$iter_warmup,
      iter_sampling = opt$iter_sampling,
      seed = opt$seed,
      adapt_delta = opt$adapt_delta,
      max_treedepth = opt$max_treedepth
    )
  )

  if (isTRUE(opt$save_full)) {
    full_path <- file.path(opt$output_dir, "full_posterior_draws.rds")
    save_object(full_object, full_path)
    message(sprintf("Saved full posterior draws to %s", full_path))
  }

  subset_registry <- vector("list", length(K_values))
  names(subset_registry) <- paste0("K", K_values)

  for (K in K_values) {
    if (K > T_total) {
      stop(sprintf("K = %d exceeds T = %d.", K, T_total), call. = FALSE)
    }

    blocks <- make_blocks(T = T_total, K = K)
    subset_fits <- vector("list", nrow(blocks))

    message(sprintf("Fitting subset posteriors for K = %d (temper = %d)...", K, K))

    for (b in seq_len(nrow(blocks))) {
      idx <- blocks$start[b]:blocks$end[b]
      block_seed <- opt$seed + 1000L * K + b

      message(sprintf(
        "  Block %d/%d: rows %d:%d (n = %d)",
        b, nrow(blocks), blocks$start[b], blocks$end[b], blocks$size[b]
      ))

      subset_fit <- fit_one_dataset(
        model = model,
        y = y[idx],
        Z = Z[idx, , drop = FALSE],
        temper = K,
        seed = block_seed,
        chains = opt$chains,
        parallel_chains = opt$parallel_chains,
        iter_warmup = opt$iter_warmup,
        iter_sampling = opt$iter_sampling,
        adapt_delta = opt$adapt_delta,
        max_treedepth = opt$max_treedepth,
        refresh = opt$refresh
      )

      subset_fits[[b]] <- list(
        block = as.list(blocks[b, , drop = FALSE]),
        seed = block_seed,
        draws = subset_fit
      )
    }

    subset_object <- list(
      K = K,
      blocks = blocks,
      fits = subset_fits,
      meta = list(
        data_file = normalizePath(opt$data, winslash = "/", mustWork = FALSE),
        stan_file = normalizePath(opt$stan_file, winslash = "/", mustWork = FALSE),
        T = T_total,
        p = p,
        temper = K,
        chains = opt$chains,
        parallel_chains = opt$parallel_chains,
        iter_warmup = opt$iter_warmup,
        iter_sampling = opt$iter_sampling,
        base_seed = opt$seed,
        adapt_delta = opt$adapt_delta,
        max_treedepth = opt$max_treedepth
      )
    )

    subset_registry[[paste0("K", K)]] <- subset_object

    if (isTRUE(opt$save_subsets)) {
      subset_path <- file.path(opt$output_dir, sprintf("subset_posteriors_K%d.rds", K))
      save_object(subset_object, subset_path)
      message(sprintf("Saved subset posterior draws to %s", subset_path))
    }
  }

  manifest <- list(
    data_file = normalizePath(opt$data, winslash = "/", mustWork = FALSE),
    stan_file = normalizePath(opt$stan_file, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE),
    K_values = K_values,
    full_file = if (isTRUE(opt$save_full)) file.path(opt$output_dir, "full_posterior_draws.rds") else NULL,
    subset_files = if (isTRUE(opt$save_subsets)) {
      setNames(
        file.path(opt$output_dir, sprintf("subset_posteriors_K%d.rds", K_values)),
        paste0("K", K_values)
      )
    } else {
      NULL
    },
    created_at = as.character(Sys.time())
  )

  manifest_path <- file.path(opt$output_dir, "fit_manifest.rds")
  save_object(manifest, manifest_path)
  message(sprintf("Saved manifest to %s", manifest_path))
}

main()
