# Fast PSD -> ACVF conversion used to avoid repeated artfimaTACVF() calls.

psd1_to_acvf <- function(P_os, df, maxlag = NULL, include_nyq = TRUE) {
  stopifnot(is.numeric(P_os), length(P_os) >= 1L, is.numeric(df), length(df) == 1L, df > 0)

  M <- length(P_os)

  if (include_nyq) {
    N <- 2 * (M - 1L)
    S0   <- P_os[1] / 2
    mid  <- if (M > 2L) P_os[2:(M - 1L)] / 2 else numeric(0)
    Snyq <- P_os[M] / 2
    S_two <- c(S0, mid, Snyq, rev(mid))
  } else {
    N <- 2L * M - 2L
    S0   <- P_os[1]
    mid  <- if (M > 1L) P_os[2:M] / 2 else numeric(0)
    S_two <- c(S0, mid, rev(mid))
  }

  acvf_full <- 2 * Re(fft(S_two, inverse = TRUE))
  acvf_full <- acvf_full * df

  out_maxlag <- if (is.null(maxlag)) {
    length(acvf_full) - 1L
  } else {
    maxlag
  }

  acvf_full[seq_len(out_maxlag + 1L)]
}

bigdog_artfima_acf <- function(n, d, lambda, phi) {
  stopifnot(
    length(n) == 1L, n >= 1,
    length(d) == 1L, length(lambda) == 1L, length(phi) == 1L
  )

  psd <- artfimaSDF(
    n = 2 * n,
    d = d,
    lambda = lambda,
    phi = phi,
    plot = "none"
  )

  psd1_to_acvf(
    P_os = psd,
    df = 1 / (2 * n),
    include_nyq = TRUE
  )[seq_len(n + 1L)]
}