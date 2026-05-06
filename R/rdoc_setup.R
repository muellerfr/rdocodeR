rdoc_default_n_nulls <- function() {
  as.integer(getOption("rdocodeR.n_nulls", 1000L))
}


rdoc_cache_root <- function() {
  tools::R_user_dir("rdocodeR", which = "cache")
}


rdoc_python_env_dir <- function() {
  file.path(rdoc_cache_root(), "python", "eigenstrapping")
}


rdoc_term_modes_file <- function() {
  path <- system.file("extdata", "rdoc_term_modes.rds", package = "rdocodeR")
  if (identical(path, "")) {
    stop("Bundled RDoC term mode table is missing from this installation.", call. = FALSE)
  }
  path
}


rdoc_term_modes <- function() {
  readRDS(rdoc_term_modes_file())
}


rdoc_eigen_surface_file <- function(hemi = c("lh", "rh")) {
  hemi <- match.arg(hemi)
  path <- system.file(
    "extdata",
    "eigenstrapping",
    paste0("fsaverage6-", hemi, ".orig.gii"),
    package = "rdocodeR"
  )
  if (identical(path, "")) {
    stop(
      sprintf("Bundled fsaverage6 surface file is missing for hemisphere '%s'.", hemi),
      call. = FALSE
    )
  }
  path
}


rdoc_split_fsavg6 <- function(x) {
  stopifnot(is.numeric(x), length(x) == 81924L)
  list(
    lh = x[seq_len(40962L)],
    rh = x[40962L + seq_len(40962L)]
  )
}


rdoc_normalize_python <- function(python) {
  if (is.null(python) || !nzchar(python)) {
    return(NULL)
  }
  normalizePath(path.expand(python), mustWork = FALSE)
}


rdoc_python_ready <- function() {
  tryCatch({
    reticulate::import("eigenstrapping", delay_load = FALSE)
    TRUE
  }, error = function(e) FALSE)
}


rdoc_use_eigen_python <- function(python = NULL, quiet = FALSE) {
  python <- rdoc_normalize_python(python)
  env_dir <- rdoc_python_env_dir()

  dir.create(dirname(env_dir), recursive = TRUE, showWarnings = FALSE)

  if (!is.null(python)) {
    if (!quiet) {
      message("Using user-specified Python interpreter: ", python)
    }
    reticulate::use_python(python, required = TRUE)
    if (!rdoc_python_ready()) {
      stop(
        paste(
          "The supplied Python interpreter does not have the 'eigenstrapping' package.",
          "Install it in that environment or rerun without `python` to let rdocodeR create a private environment."
        ),
        call. = FALSE
      )
    }
    return(invisible(list(method = "python", path = python)))
  }

  if (reticulate::virtualenv_exists(env_dir)) {
    if (!quiet) {
      message("Using cached Python environment: ", env_dir)
    }
    reticulate::use_virtualenv(env_dir, required = TRUE)
  } else {
    starter <- Sys.which("python3")
    if (!nzchar(starter)) {
      stop(
        paste(
          "No system 'python3' interpreter was found to create the private eigenstrapping environment.",
          "Install Python 3 or call `rdoc_setup(python = \"/path/to/python\")` with an interpreter that has eigenstrapping available."
        ),
        call. = FALSE
      )
    }
    if (!quiet) {
      message("Creating private Python environment at: ", env_dir)
    }
    reticulate::virtualenv_create(
      envname = env_dir,
      python = starter,
      packages = NULL
    )
    reticulate::use_virtualenv(env_dir, required = TRUE)
  }

  if (!rdoc_python_ready()) {
    if (!quiet) {
      message("Installing Python package 'eigenstrapping' into: ", env_dir)
    }
    reticulate::py_install(
      packages = "eigenstrapping",
      envname = env_dir,
      method = "virtualenv",
      pip = TRUE
    )
  }

  if (!rdoc_python_ready()) {
    stop(
      paste(
        "The private Python environment was created, but 'eigenstrapping' could not be imported.",
        "Please rerun with `python = \"/path/to/python\"` if you already have a working environment."
      ),
      call. = FALSE
    )
  }

  invisible(list(method = "virtualenv", path = env_dir))
}


rdoc_create_eigenstraps <- function(fs_overlay, n_nulls, num_modes, n_jobs = 1L) {
  hemis <- c("lh", "rh")
  np <- reticulate::import("numpy")
  es <- reticulate::import("eigenstrapping")

  Sys.setenv(MPLCONFIGDIR = file.path(tempdir(), "matplotlib"))

  surface_files <- setNames(vapply(
    hemis,
    rdoc_eigen_surface_file,
    character(1)
  ), hemis)

  surf_eigen <- setNames(lapply(seq_along(hemis), function(i) {
    hemi <- hemis[[i]]
    es$SurfaceEigenstrapping(
      data = np$array(fs_overlay[[hemi]]),
      surface = surface_files[[hemi]],
      num_modes = as.integer(num_modes[[i]]),
      resample = TRUE,
      n_jobs = as.integer(n_jobs)
    )
  }), hemis)

  setNames(lapply(hemis, function(hemi) {
    surf_eigen[[hemi]](n = as.integer(n_nulls))
  }), hemis)
}


rdoc_expected_term_null_files <- function(term_ids, n_nulls = rdoc_default_n_nulls()) {
  dir <- rdoc_term_nulls_dir(n_nulls = n_nulls, create = FALSE)
  file.path(dir, paste0(term_ids, "_eigenstraps.rds"))
}


#' Local Cache Directory for Generated RDoC Term Null Maps
#'
#' Returns the user cache directory that stores generated per-term null map files.
#'
#' @param n_nulls Number of cached nulls per term.
#' @param create Logical; if `TRUE`, create the directory if needed.
#'
#' @return Absolute file path to the cache directory used for the requested null-bank size.
#' @export
rdoc_term_nulls_dir <- function(n_nulls = rdoc_default_n_nulls(), create = FALSE) {
  n_nulls <- as.integer(n_nulls[[1]])
  if (!is.finite(n_nulls) || n_nulls < 1L) {
    stop("`n_nulls` must be a positive integer.", call. = FALSE)
  }

  path <- file.path(rdoc_cache_root(), "term_nulls", paste0("n", n_nulls))
  if (isTRUE(create)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}


#' Set Up Cached RDoC Term Null Maps
#'
#' Creates or refreshes the cached per-term eigenstrap null maps used by [rdoc_decode()].
#' These null maps are generated locally for the user and are not shipped precomputed with the
#' package. By default, the setup generates `1000` nulls per bundled RDoC term map and stores
#' them in the user cache directory returned by [rdoc_term_nulls_dir()].
#'
#' @param n_nulls Number of eigenstraps to generate per term. Default is `1000`.
#' @param force Logical; if `TRUE`, regenerate files even when they already exist.
#' @param python Optional Python interpreter path. If `NULL`, `rdocodeR` creates and reuses
#'   a private virtual environment in the user cache directory.
#' @param terms Optional character vector of term ids to generate. Defaults to all bundled terms.
#' @param n_jobs Number of Python jobs used per hemisphere during eigenstrap generation.
#' @param quiet Logical; suppress progress messages.
#'
#' @return Invisibly returns a list with cache metadata and generated files.
#' @export
rdoc_setup <- function(n_nulls = rdoc_default_n_nulls(),
                       force = FALSE,
                       python = NULL,
                       terms = NULL,
                       n_jobs = 1L,
                       quiet = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      paste(
        "The `reticulate` package is required to generate cached term null maps.",
        "Please reinstall `rdocodeR` with its dependencies."
      ),
      call. = FALSE
    )
  }

  n_nulls <- as.integer(n_nulls[[1]])
  n_jobs <- as.integer(n_jobs[[1]])
  if (!is.finite(n_nulls) || n_nulls < 1L) {
    stop("`n_nulls` must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(n_jobs) || n_jobs < 1L) {
    stop("`n_jobs` must be a positive integer.", call. = FALSE)
  }

  terms_list <- readRDS(rdoc_terms_file())
  term_ids <- names(terms_list)

  if (is.null(term_ids) || any(!nzchar(term_ids))) {
    stop("Bundled term maps must be a named list.", call. = FALSE)
  }

  if (!is.null(terms)) {
    terms <- unique(as.character(terms))
    missing_terms <- setdiff(terms, term_ids)
    if (length(missing_terms)) {
      stop(
        "Requested terms are not available in the bundled RDoC term maps: ",
        paste(missing_terms, collapse = ", "),
        call. = FALSE
      )
    }
    terms_list <- terms_list[terms]
    term_ids <- names(terms_list)
  }

  mode_table <- rdoc_term_modes()
  missing_mode_rows <- setdiff(term_ids, mode_table$term)
  if (length(missing_mode_rows)) {
    stop(
      "Bundled term mode table is missing terms: ",
      paste(missing_mode_rows, collapse = ", "),
      call. = FALSE
    )
  }
  mode_table <- mode_table[match(term_ids, mode_table$term), , drop = FALSE]

  cache_dir <- rdoc_term_nulls_dir(n_nulls = n_nulls, create = TRUE)

  files <- file.path(cache_dir, paste0(term_ids, "_eigenstraps.rds"))
  needs_generation <- force | !file.exists(files)

  if (!any(needs_generation)) {
    if (!quiet) {
      message("Cached term null maps already available at: ", cache_dir)
    }
    return(invisible(list(
      n_nulls = n_nulls,
      cache_dir = cache_dir,
      generated = character(0),
      reused = basename(files)
    )))
  }

  rdoc_use_eigen_python(python = python, quiet = quiet)

  generated <- character(0)
  reused <- character(0)

  for (i in seq_along(term_ids)) {
    term_id <- term_ids[[i]]
    target_file <- files[[i]]

    if (!needs_generation[[i]]) {
      reused <- c(reused, basename(target_file))
      next
    }

    modes <- c(mode_table$best_mode_lh[[i]], mode_table$best_mode_rh[[i]])
    if (!quiet) {
      message(
        sprintf(
          "[%d/%d] %s | modes: lh=%d rh=%d | nulls=%d",
          i,
          length(term_ids),
          term_id,
          modes[[1]],
          modes[[2]],
          n_nulls
        )
      )
    }

    term_nulls <- rdoc_create_eigenstraps(
      fs_overlay = rdoc_split_fsavg6(terms_list[[i]]),
      n_nulls = n_nulls,
      num_modes = modes,
      n_jobs = n_jobs
    )

    tmp_file <- tempfile(
      pattern = paste0(term_id, "_"),
      tmpdir = cache_dir,
      fileext = ".rds"
    )
    saveRDS(term_nulls, tmp_file, compress = "xz")
    ok <- file.rename(tmp_file, target_file)
    if (!isTRUE(ok)) {
      stop("Failed to move generated null map file into the cache directory.", call. = FALSE)
    }

    generated <- c(generated, basename(target_file))
    gc()
  }

  invisible(list(
    n_nulls = n_nulls,
    cache_dir = cache_dir,
    generated = generated,
    reused = reused
  ))
}
