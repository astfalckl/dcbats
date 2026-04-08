
# ------------------------------------------------------------------------------
# This script reproduces the values in Table 3
# ------------------------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(data.table)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

fit_chunk <- function(
    Y_chunk, K, sigma1_0 = 1, sigma2_0 = 1, chains = 4, iter_warmup = 5000,
    iter_sampling = 5000, seed = 123
){
  data_list <- list(
    T = nrow(Y_chunk),
    dim = ncol(Y_chunk),
    Y = lapply(seq_len(nrow(Y_chunk)), function(t) Y_chunk[t, ]),
    power = K,
    sigma1_0 = sigma1_0,
    sigma2_0 = sigma2_0
  )

  mod$sample(
    data = data_list,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = 200
  )
}

get_draws <- function(fit, vars) {
  as_draws_df(fit$draws(variables = vars))
}


# ------------------------------------------------------------------------------
# Fit Stan Models
# ------------------------------------------------------------------------------

Y <- as.matrix(fread("la_pm_log_positive.csv"))
T <- nrow(Y)
dim <- ncol(Y)

K <- 10
n_chunk <- T %/% K

chunks <- lapply(seq_len(K), function(k) {
  idx <- ((k - 1) * n_chunk + 1):(k * n_chunk)
  Y[idx, , drop = FALSE]
})

mod <- cmdstan_model("dcc_garch.stan")

chunk_fits <- lapply(seq_along(chunks), function(k) {
  fit_chunk(chunks[[k]], K = K, seed = 100 + k)
})

full_fit <- mod$sample(
  data = list(
    T = nrow(Y),
    dim = ncol(Y),
    Y = lapply(seq_len(nrow(Y)), function(t) Y[t, ]),
    power = 1,
    sigma1_0 = 1,
    sigma2_0 = 1
  ),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  seed = 999
)

# ------------------------------------------------------------------------------
# Aggregate Results
# ------------------------------------------------------------------------------

params <- c("a1", "a2", "b1", "b2", "w1", "w2", "mu1", "mu2", "r")

chunk_draws <- lapply(chunk_fits, get_draws, vars = params)
full_draws  <- get_draws(full_fit, vars = params)

ci_full <- rbindlist(lapply(params, function(p) {
  qs <- quantile(full_draws[[p]], probs = c(0.025, 0.975), names = FALSE)
  data.table(parameter = p, lower = qs[1], upper = qs[2])
}))

ci_dc <- rbindlist(lapply(params, function(p) {
    qlo <- mean(vapply(chunk_draws, function(d) quantile(d[[p]], 0.025, names = FALSE), numeric(1)))
    qhi <- mean(vapply(chunk_draws, function(d) quantile(d[[p]], 0.975, names = FALSE), numeric(1)))
    data.table(parameter = p, lower = qlo, upper = qhi)
}))

result <- merge(ci_dc, ci_full, by = "parameter", suffixes = c("_dc", "_full"))
print(result)





















get_draws <- function(fit, vars) {
  as_draws_df(fit$draws(variables = vars))
}

chunk_draws <- lapply(chunk_fits, get_draws, vars = params)
full_draws  <- get_draws(full_fit, vars = params)

wasserstein_barycenter_1d <- function(draw_list, probs = seq(0.001, 0.999, by = 0.001)) {
  out <- vector("list", length = length(params))
  names(out) <- params

  for (p in params) {
    qmat <- sapply(draw_list, function(d) quantile(d[[p]], probs = probs, names = FALSE))
    qbar <- rowMeans(qmat)
    out[[p]] <- data.frame(prob = probs, qbar = qbar)
  }
  out
}

wbary <- wasserstein_barycenter_1d(chunk_draws)


sample_barycenter_draws <- function(wbary_obj, M = 4000) {
  u <- sort(runif(M))
  out <- data.frame(.draw = seq_len(M))
  for (p in names(wbary_obj)) {
    out[[p]] <- approx(
      x = wbary_obj[[p]]$prob,
      y = wbary_obj[[p]]$qbar,
      xout = u,
      rule = 2
    )$y
  }
  out
}

bary_draws <- sample_barycenter_draws(wbary, M = 4000)


ci_from_draws <- function(draw_df, vars, level = 0.95) {
  alpha <- (1 - level) / 2
  rbindlist(lapply(vars, function(v) {
    qs <- quantile(draw_df[[v]], probs = c(alpha, 1 - alpha), names = FALSE)
    data.table(parameter = v, lower = qs[1], upper = qs[2])
  }))
}

ci_full <- ci_from_draws(full_draws, params)
ci_dc   <- ci_from_draws(bary_draws, params)

result <- merge(ci_dc, ci_full, by = "parameter", suffixes = c("_dc", "_full"))
print(result)
