rdoc_spatial_cor_with_nulls_pairwise <- function(x,
                                                 y,
                                                 x_nulls,
                                                 method = "pearson",
                                                 use = "pairwise.complete.obs") {
  x_nulls <- rdoc_normalize_null_bank(x_nulls, n_vertices = length(y))

  n_nulls <- nrow(x_nulls)
  r_obs <- stats::cor(x, y, use = use, method = method)
  r_null <- apply(
    x_nulls,
    1L,
    function(x_null) stats::cor(x_null, y, use = use, method = method)
  )

  if (!is.finite(r_obs)) {
    p_perm <- NA_real_
  } else {
    p_perm <- (sum(abs(r_null) >= abs(r_obs), na.rm = TRUE) + 1) / (1 + n_nulls)
  }

  data.frame(r.obs = r_obs, p.perm = p_perm, stringsAsFactors = FALSE)
}


rdoc_normalize_null_bank <- function(x_nulls, n_vertices) {
  if (is.list(x_nulls) && all(c("lh", "rh") %in% names(x_nulls))) {
    if (is.null(dim(x_nulls$lh)) || length(dim(x_nulls$lh)) == 1L) {
      x_nulls$lh <- matrix(x_nulls$lh, nrow = 1L)
    }
    if (is.null(dim(x_nulls$rh)) || length(dim(x_nulls$rh)) == 1L) {
      x_nulls$rh <- matrix(x_nulls$rh, nrow = 1L)
    }
    x_nulls <- cbind(x_nulls$lh, x_nulls$rh)
  }

  if (!is.matrix(x_nulls)) {
    x_nulls <- as.matrix(x_nulls)
  }

  if (ncol(x_nulls) == n_vertices) {
    return(x_nulls)
  }
  if (nrow(x_nulls) == n_vertices) {
    return(t(x_nulls))
  }

  stop(
    sprintf(
      "Null bank must have %d vertices in either rows or columns.",
      as.integer(n_vertices)
    ),
    call. = FALSE
  )
}


rdoc_term_null_file <- function(term_id, n_nulls = rdoc_default_n_nulls()) {
  path <- file.path(rdoc_term_nulls_dir(n_nulls = n_nulls, create = FALSE), paste0(term_id, "_eigenstraps.rds"))
  if (!file.exists(path)) {
    stop(
      sprintf("Cached term null file is missing for term '%s'.", term_id),
      call. = FALSE
    )
  }
  path
}


#' Decode a Brain Overlay Against RDoC Term Maps
#'
#' Computes spatial correlations between a cortical overlay and the internal RDoC term maps,
#' using cached per-term eigenstrap null maps to estimate permutation p-values. These null maps
#' are not shipped precomputed with the package. By default, `rdocodeR` will generate the
#' required null maps on first use and then reuse them from the user cache directory on
#' subsequent runs.
#'
#' @param fs_overlay Numeric cortical overlay vector matching the bundled RDoC term maps.
#' @param save_results Logical; write a TSV with correlations.
#' @param results_file Optional output path for the TSV file when `save_results = TRUE`.
#'   If `NULL`, defaults to `rdoc_decode_results.tsv` in the current working directory.
#' @param cor_method Correlation method (`"pearson"` or `"spearman"`).
#'   Default is `"pearson"`.
#' @param absolute_r Logical; if `TRUE`, stores absolute values of `r`.
#' @param n_nulls Number of cached eigenstrap nulls per term. Default is `1000`.
#' @param auto_setup Logical; if `TRUE`, run [rdoc_setup()] automatically when cached nulls
#'   are missing.
#' @param python Optional Python interpreter path forwarded to [rdoc_setup()] when
#'   `auto_setup = TRUE`.
#' @param setup_n_jobs Number of Python jobs used during automatic null generation.
#' @param mc_cores Number of cores used for Unix `mclapply`.
#'
#' @return A data frame with columns `Domain`, `Term`, `r`, and `p`.
#' @export
rdoc_decode <- function(fs_overlay,
                        save_results = FALSE,
                        results_file = NULL,
                        cor_method = c("pearson", "spearman"),
                        absolute_r = FALSE,
                        n_nulls = rdoc_default_n_nulls(),
                        auto_setup = TRUE,
                        python = NULL,
                        setup_n_jobs = 1L,
                        mc_cores = max(1L, parallel::detectCores() - 1L)) {
  cor_method <- match.arg(cor_method)
  n_nulls <- as.integer(n_nulls[[1]])
  if (!is.finite(n_nulls) || n_nulls < 1L) {
    stop("`n_nulls` must be a positive integer.", call. = FALSE)
  }

  terms <- readRDS(rdoc_terms_file())
  ref <- rdoc_terms_reference()

  if (!is.numeric(fs_overlay) || length(dim(fs_overlay)) > 1L) {
    stop("`fs_overlay` must be a numeric vector.", call. = FALSE)
  }
  if (length(terms) != nrow(ref)) {
    stop(
      sprintf(
        "Bundled term map file contains %d maps, expected %d.",
        length(terms),
        nrow(ref)
      ),
      call. = FALSE
    )
  }
  if (length(fs_overlay) != length(terms[[1L]])) {
    stop(
      sprintf(
        "`fs_overlay` has length %d, expected %d vertices.",
        length(fs_overlay),
        length(terms[[1L]])
      ),
      call. = FALSE
    )
  }

  term_ids <- names(terms)
  if (is.null(term_ids) || any(!nzchar(term_ids))) {
    term_ids <- sprintf("term_%03d", seq_along(terms))
  }

  null_files <- rdoc_expected_term_null_files(term_ids, n_nulls = n_nulls)
  if (isTRUE(auto_setup) && any(!file.exists(null_files))) {
    setup_res <- rdoc_setup(
      n_nulls = n_nulls,
      python = python,
      n_jobs = setup_n_jobs
    )
    null_files <- file.path(setup_res$cache_dir, paste0(term_ids, "_eigenstraps.rds"))
  }
  if (any(!file.exists(null_files))) {
    stop(
      paste(
        "Cached RDoC term null maps are missing.",
        "Run `rdoc_setup()` first or call `rdoc_decode(..., auto_setup = TRUE)`."
      ),
      call. = FALSE
    )
  }

  worker <- function(i) {
    test <- tryCatch(
      rdoc_spatial_cor_with_nulls_pairwise(
        x = fs_overlay,
        y = terms[[i]],
        x_nulls = readRDS(null_files[[i]]),
        method = cor_method
      ),
      error = function(e) NULL
    )

    if (is.null(test)) {
      return(data.frame(
        Domain = ref$Domain[i],
        Term = ref$Term[i],
        r = NA_real_,
        p = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      Domain = ref$Domain[i],
      Term = ref$Term[i],
      r = as.numeric(test$r.obs),
      p = as.numeric(test$p.perm),
      stringsAsFactors = FALSE
    )
  }

  idx <- seq_along(terms)
  if (.Platform$OS.type == "windows") {
    rows <- lapply(idx, worker)
  } else {
    rows <- parallel::mclapply(idx, worker, mc.cores = as.integer(mc_cores))
  }

  corr_df <- do.call(rbind, rows)

  if (absolute_r) {
    corr_df$r <- abs(corr_df$r)
  }

  if (isTRUE(save_results)) {
    if (is.null(results_file) || identical(results_file, "")) {
      results_file <- file.path(getwd(), "rdoc_decode_results.tsv")
    }
    dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      corr_df,
      file = results_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  corr_df
}
