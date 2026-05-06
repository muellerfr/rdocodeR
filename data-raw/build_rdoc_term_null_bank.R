#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(reticulate)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

split_fsavg6 <- function(x) {
  stopifnot(is.numeric(x), length(x) == 81924L)
  list(
    lh = x[seq_len(40962L)],
    rh = x[40962L + seq_len(40962L)]
  )
}

create_eigenstraps_fixed <- function(fs.overlay, n.nulls, num.modes) {
  hemis <- c("lh", "rh")
  np <- reticulate::import("numpy")
  es <- reticulate::import("eigenstrapping")

  Sys.setenv(MPLCONFIGDIR = file.path(tempdir(), "matplotlib"))

  surface.filename <- setNames(lapply(hemis, function(h) {
    file.path(repo_root, "inst", "extdata", "eigenstrapping", paste0("fsaverage6-", h, ".orig.gii"))
  }), hemis)

  surf_eigen <- setNames(lapply(seq_along(hemis), function(i) {
    h <- hemis[[i]]
    es$SurfaceEigenstrapping(
      data = np$array(fs.overlay[[h]]),
      surface = surface.filename[[h]],
      num_modes = as.integer(num.modes[[i]]),
      resample = TRUE,
      n_jobs = 10L
    )
  }), hemis)

  setNames(lapply(hemis, function(h) {
    surf_eigen[[h]](n = as.integer(n.nulls))
  }), hemis)
}

normalize_summary <- function(x) {
  if (is.data.frame(x) && all(c("term", "best_mode_lh", "best_mode_rh") %in% names(x))) {
    return(x[, c("term", "best_mode_lh", "best_mode_rh")])
  }
  data.frame(
    term = names(x),
    best_mode_lh = vapply(x, function(z) z$best_mode_lh, numeric(1)),
    best_mode_rh = vapply(x, function(z) z$best_mode_rh, numeric(1)),
    stringsAsFactors = FALSE
  )
}

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- if (length(script_file)) {
  normalizePath(sub("^--file=", "", script_file[[1]]))
} else {
  file.path(getwd(), "data-raw", "build_rdoc_term_null_bank.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) {
    return(default)
  }
  sub(paste0("^", flag, "="), "", hit[[1]])
}

terms_file <- arg_value(
  "--terms-file",
  file.path(repo_root, "inst", "extdata", "rdoc_terms.fsavg6.lh.rh.rds")
)
summary_file <- arg_value(
  "--summary-file",
  file.path(repo_root, "inst", "extdata", "rdoc_term_modes.rds")
)
per_term_dir <- arg_value(
  "--per-term-dir",
  file.path(repo_root, "data-raw", "rdoc_term_nulls_2000")
)
n_nulls <- as.integer(arg_value("--n-nulls", "2000"))
selected_terms <- arg_value("--terms", NULL)
python_path <- arg_value("--python", NULL)
if (!is.null(selected_terms) && nzchar(selected_terms)) {
  selected_terms <- strsplit(selected_terms, ",", fixed = TRUE)[[1]]
  selected_terms <- trimws(selected_terms)
  selected_terms <- selected_terms[nzchar(selected_terms)]
}
if (!is.null(python_path) && nzchar(python_path)) {
  python_path <- path.expand(python_path)
  if (grepl("/bin/python([0-9.]*)?$", python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  } else {
    reticulate::use_condaenv(python_path, required = TRUE)
  }
}

if (!file.exists(terms_file)) {
  stop("Term map file not found: ", terms_file, call. = FALSE)
}
if (!file.exists(summary_file)) {
  stop("Summary file not found: ", summary_file, call. = FALSE)
}

dir.create(per_term_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading packaged term maps...")
terms <- readRDS(terms_file)
summary_df <- normalize_summary(readRDS(summary_file))

if (length(terms) != nrow(summary_df)) {
  stop(
    sprintf(
      "Term map count (%d) does not match summary count (%d).",
      length(terms),
      nrow(summary_df)
    ),
    call. = FALSE
  )
}

missing_terms <- setdiff(names(terms), summary_df$term)
if (length(missing_terms)) {
  stop(
    "Missing terms in summary: ",
    paste(missing_terms, collapse = ", "),
    call. = FALSE
  )
}

summary_df <- summary_df[match(names(terms), summary_df$term), ]

if (!is.null(selected_terms)) {
  missing_selected <- setdiff(selected_terms, names(terms))
  if (length(missing_selected)) {
    stop(
      "Requested terms not found in packaged data: ",
      paste(missing_selected, collapse = ", "),
      call. = FALSE
    )
  }
  terms <- terms[selected_terms]
  summary_df <- summary_df[match(names(terms), summary_df$term), ]
}

for (i in seq_along(terms)) {
  term_name <- names(terms)[[i]]
  term_map <- terms[[i]]
  modes <- c(summary_df$best_mode_lh[[i]], summary_df$best_mode_rh[[i]])

  message(
    sprintf(
      "[%d/%d] %s | modes: lh=%d rh=%d | nulls=%d",
      i,
      length(terms),
      term_name,
      modes[[1]],
      modes[[2]],
      n_nulls
    )
  )

  term_nulls <- create_eigenstraps_fixed(
    fs.overlay = split_fsavg6(term_map),
    n.nulls = n_nulls,
    num.modes = modes
  )

  per_term_file <- file.path(per_term_dir, paste0(term_name, "_eigenstraps.rds"))
  saveRDS(term_nulls, per_term_file, compress = "xz")

  gc()
}

message("Done.")
message("Per-term nulls: ", per_term_dir)
