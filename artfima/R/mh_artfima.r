# Custom Metropolis-Hastings sampler for ARTFIMA model.
# Kept close to the working research implementation.

mh_artfima <- function(
  x,
  n_iter = 1000,
  init = list(d = 0.4, lambda = 0.01, phi = 0.99, sigma2 = 1),
  sd_prop = list(d = 0.01, lambda = 0.001, phi = 0.001, sigma2 = 0.01),
  incl_sample = list(
    sample_d = TRUE,
    sample_lambda = TRUE,
    sample_phi = TRUE,
    sample_sigma2 = TRUE
  ),
  temper = 1
) {
  acf_n <- length(x) - 1L

  logpost <- function(par) {
    if (par$d <= 0 ||
        par$d >= 0.5 ||
        par$lambda <= 0 ||
        abs(par$phi) >= 1 ||
        par$sigma2 <= 0) {
      return(-Inf)
    }

    acf <- bigdog_artfima_acf(
      n = acf_n,
      d = par$d,
      lambda = par$lambda,
      phi = par$phi
    )

    temper * dnormtz(
      x,
      mu = 0,
      acf = acf,
      log = TRUE,
      method = "gschur"
    )
  }

  out <- matrix(
    NA_real_,
    nrow = n_iter,
    ncol = 4L,
    dimnames = list(NULL, c("d", "lambda", "phi", "sigma2"))
  )

  cur <- init
  lcur <- logpost(cur)
  acc <- 0L

  prop_d <- cur$d
  prop_lambda <- cur$lambda
  prop_phi <- cur$phi
  prop_sigma2 <- cur$sigma2

  for (i in seq_len(n_iter)) {
    cat(sprintf("\rIteration %d/%d", i, n_iter))

    if (isTRUE(incl_sample$sample_d)) {
      prop_d <- cur$d + rnorm(1, 0, sd_prop$d)
    }

    if (isTRUE(incl_sample$sample_lambda)) {
      prop_lambda <- 10^(log10(cur$lambda) + rnorm(1, 0, sd_prop$lambda))
    }

    if (isTRUE(incl_sample$sample_phi)) {
      prop_phi <- cur$phi + rnorm(1, 0, sd_prop$phi)
    }

    if (isTRUE(incl_sample$sample_sigma2)) {
      prop_sigma2 <- cur$sigma2 + rnorm(1, 0, sd_prop$sigma2)
    }

    prop <- list(
      d = prop_d,
      lambda = prop_lambda,
      phi = prop_phi,
      sigma2 = prop_sigma2
    )

    lprop <- logpost(prop)

    if (log(runif(1)) < (lprop - lcur)) {
      cur <- prop
      lcur <- lprop
      acc <- acc + 1L
    }

    out[i, ] <- c(cur$d, cur$lambda, cur$phi, cur$sigma2)
  }

  cat("\n")

  attr(out, "accept_rate") <- acc / n_iter
  out
}