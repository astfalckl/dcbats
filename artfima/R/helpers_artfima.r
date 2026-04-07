# Shared package setup and sourcing for ARTFIMA code.

required_pkgs <- c(
  "ltsa",
  "devtools",
  "gsl",
  "SuperGauss",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "readr",
  "tibble"
)

load_required_packages <- function(pkgs = required_pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Missing required packages: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(lapply(pkgs, library, character.only = TRUE))
}

source_artfima_code <- function(artfima_dir = "artfima/R") {
  if (!dir.exists(artfima_dir)) {
    stop(sprintf("Directory not found: %s", artfima_dir), call. = FALSE)
  }

  artfima_files <- list.files(
    artfima_dir,
    full.names = TRUE,
    pattern = "\\.[Rr]$"
  )

  if (length(artfima_files) == 0) {
    stop(sprintf("No R files found in: %s", artfima_dir), call. = FALSE)
  }

  invisible(lapply(artfima_files, source))
}

ensure_project_subdirs <- function(paths) {
  invisible(vapply(paths, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

make_scenario_id <- function(lambda, d, phi) {
  paste0(
    "lambda", gsub("\\.", "", format(lambda, scientific = FALSE, trim = TRUE)),
    "_d", gsub("\\.", "", format(d, scientific = FALSE, trim = TRUE)),
    "_phi", gsub("\\.", "", format(phi, scientific = FALSE, trim = TRUE))
  )
}

make_sim_id <- function(i, width = 3) {
  stringr::str_pad(i, width = width, side = "left", pad = "0")
}

default_param_grid <- function() {
  expand.grid(
    lambda = c(0.005, 0.1),
    phi    = c(0.1, 0.9, 0.99),
    d      = c(0.1, 0.3)
  )
}

default_mcmc_proposals <- function(phi, d) {
  phi_prop <- dplyr::case_when(
    phi == 0.1  ~ 0.02,
    phi == 0.9  ~ 0.008,
    phi == 0.99 ~ 0.002,
    TRUE ~ NA_real_
  )

  d_prop <- dplyr::case_when(
    d == 0.1 ~ 0.02,
    d == 0.3 ~ 0.02,
    TRUE ~ NA_real_
  )

  if (is.na(phi_prop) || is.na(d_prop)) {
    stop("No default proposal scale defined for supplied phi/d.", call. = FALSE)
  }

  list(
    d = d_prop,
    lambda = 0.3,
    phi = phi_prop,
    sigma2 = 0.01
  )
}

save_metadata_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
}