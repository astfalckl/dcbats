compute_percentiles <- function(x, probs = seq(0.01, 0.99, by = 0.01)) {
  quantile(x, probs = probs, type = 7, names = FALSE)
}

wasserstein_from_percentiles <- function(
  qP, qQ, p = 2, method = c("trap", "riemann")
) {
  stopifnot(length(qP) == length(qQ), p >= 1)

  method <- match.arg(method)
  du <- 1 / (length(qP) + 1)
  d <- abs(qP - qQ)^p

  wp_p <- if (method == "trap") {
    du * (0.5 * d[1] + sum(d[2:(length(d) - 1)]) + 0.5 * d[length(d)])
  } else {
    du * sum(d)
  }

  wp_p^(1 / p)
}

W1_rel <- function(qP, qQ, method = c("trap", "riemann")) {
  method <- match.arg(method)

  num <- wasserstein_from_percentiles(qP, qQ, p = 1, method = method)

  mid <- ceiling(length(qQ) / 2)
  m <- qQ[mid]
  a <- abs(qQ - m)
  du <- 1 / (length(qQ) + 1)

  denom <- if (method == "trap") {
    du * (0.5 * a[1] + sum(a[2:(length(a) - 1)]) + 0.5 * a[length(a)])
  } else {
    du * sum(a)
  }

  num / denom
}

average_subset_quantiles <- function(df, value_col = "value") {
  value_sym <- rlang::sym(value_col)

  df %>%
    dplyr::group_by(.data$p) %>%
    dplyr::summarise(value = mean(!!value_sym), .groups = "drop")
}