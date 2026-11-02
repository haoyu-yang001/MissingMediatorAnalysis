#' One-row summary of a mediation fit
#'
#' @param x an object of class \code{"mma"}.
#' @param row.names,optional ignored.
#' @param ... ignored.
#' @return A one-row data frame with sample sizes, path estimates, effects,
#'   tests and, for sequential SynSurr, the decomposition of
#'   \eqn{\hat\beta_M} into its direct and surrogate parts and the partial
#'   \eqn{R^2} of the mediator on the surrogates.
#' @export
as.data.frame.mma <- function(x, row.names = NULL, optional = FALSE, ...) {
  e <- x$effects
  g <- function(r, c) e[r, c]
  st <- x$fits$outcome$steps
  data.frame(
    exposure = x$vars$exposure, mediator = x$vars$mediator,
    method_mediator = x$methods[["mediator"]], method_outcome = x$methods[["outcome"]],
    n = x$n, n_M = x$n_obs, n_Y = x$n_obs_Y, n_both = x$n_obs_both,
    alpha = g("alpha", "estimate"), se_alpha = g("alpha", "se"), p_alpha = g("alpha", "p"),
    beta = g("beta", "estimate"), se_beta = g("beta", "se"), p_beta = g("beta", "p"),
    nie = g("NIE", "estimate"), se_nie = g("NIE", "se"),
    p_sobel = x$tests$p[1], p_joint = x$tests$p[2], p_mc = x$tests$p[3],
    mc_lo = if (nrow(x$nie_ci) > 1) x$nie_ci$lower[2] else NA,
    mc_hi = if (nrow(x$nie_ci) > 1) x$nie_ci$upper[2] else NA,
    nde = g("NDE", "estimate"), se_nde = g("NDE", "se"), p_nde = g("NDE", "p"),
    te = g("TE", "estimate"), se_te = g("TE", "se"), p_te = g("TE", "p"),
    pm = g("PM", "estimate"), cor_ab = x$cor_ab,
    beta_direct = if (is.null(st)) NA else st$beta_M_direct,
    beta_surrogate = if (is.null(st)) NA else st$beta_M_surrogate,
    partial_r2 = if (is.null(st)) NA else st$partial_r2,
    row.names = NULL)
}

#' Mediation analysis over many exposures and mediators
#'
#' Calls \code{\link{mma}} for every exposure x mediator pair and returns one
#' row per pair. Failures and warnings are recorded in the columns
#' \code{status} and \code{message} instead of stopping the scan. Rows with a
#' missing exposure are dropped pair by pair.
#'
#' @param data a data frame with the mediators, outcome, surrogates and
#'   covariates.
#' @param exposures either a character vector of exposure columns of
#'   \code{data}, or a numeric matrix with \code{nrow(data)} rows and column
#'   names (e.g. thousands of proteins), which avoids copying them into
#'   \code{data}.
#' @param mediators character vector of mediator columns.
#' @param surrogates character vector of the surrogate columns, one per
#'   mediator (in the same order).
#' @param outcome,covariates,outcome_surrogate as in \code{\link{mma}}.
#' @param auxiliary optional fully observed auxiliary columns added to every
#'   mediator's surrogate (sequential SynSurr only; see \code{\link{mma}}).
#' @param ... further arguments to \code{\link{mma}} (e.g.
#'   \code{outcome_method}, \code{mediator_method}, \code{mc_draws}).
#' @param n_cores number of cores (forked with \code{parallel::mclapply}; 1
#'   runs serially). Exposures whose worker fails (e.g. out of memory) are
#'   returned as error rows.
#' @param seed optional integer; if given, the random numbers of the Monte Carlo
#'   interval for the \eqn{j}-th exposure use seed \code{seed + j}, so results
#'   are reproducible also in parallel.
#' @return A data frame with columns \code{exposure}, \code{mediator},
#'   \code{status} (\code{"ok"}, \code{"warning"} or \code{"error"}),
#'   \code{message}, and the columns of \code{\link{as.data.frame.mma}}.
#' @examples
#' dat <- sim_mma_data(n = 2000, miss_rate = 0.75, rho = 0.8,
#'                     rho_Y = 0.8, miss_rate_Y = 0.5)
#' G <- cbind(A = dat$A, A2 = rnorm(2000))
#' mma_scan(dat, G, mediators = "M", surrogates = "S", outcome = "Y",
#'          covariates = "X", outcome_surrogate = "S_Y", mc_draws = 0)
#' @export
mma_scan <- function(data, exposures, mediators, surrogates, outcome,
                     covariates = NULL, outcome_surrogate = NULL, auxiliary = NULL, ...,
                     n_cores = 1, seed = NULL) {
  data <- as.data.frame(data)
  stopifnot(length(surrogates) == length(mediators))
  if (is.character(exposures)) {
    G <- as.matrix(data[, exposures, drop = FALSE])
  } else {
    G <- as.matrix(exposures)
    if (nrow(G) != nrow(data) || is.null(colnames(G)))
      stop("an exposure matrix needs nrow(data) rows and column names")
  }
  one <- function(j) {
    if (!is.null(seed)) set.seed(seed + j)
    keep <- !is.na(G[, j])
    d <- data[keep, , drop = FALSE]
    d$.exposure <- G[keep, j]
    do.call(rbind, lapply(seq_along(mediators), function(i) {
      msg <- character()
      fit <- withCallingHandlers(
        tryCatch(mma(d, exposure = ".exposure", mediator = mediators[i], outcome = outcome,
                     surrogate = c(surrogates[i], auxiliary), covariates = covariates,
                     outcome_surrogate = outcome_surrogate, keep_IF = FALSE,
                     verbose = FALSE, ...),
                 error = function(e) e),
        warning = function(w) {
          msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning")
        })
      head <- data.frame(exposure = colnames(G)[j], mediator = mediators[i])
      if (inherits(fit, "error"))
        return(cbind(head, status = "error", message = conditionMessage(fit)))
      row <- as.data.frame(fit)
      row$exposure <- NULL; row$mediator <- NULL
      cbind(head, status = if (length(msg)) "warning" else "ok",
            message = paste(msg, collapse = "; "), row)
    }))
  }
  rows <- if (n_cores > 1) parallel::mclapply(seq_len(ncol(G)), one, mc.cores = n_cores) else
    lapply(seq_len(ncol(G)), one)
  failed <- which(!vapply(rows, is.data.frame, logical(1)))
  for (j in failed)
    rows[[j]] <- data.frame(exposure = colnames(G)[j], mediator = mediators, status = "error",
                            message = paste("worker failed:", paste(as.character(rows[[j]]), collapse = " ")))
  .rbind_fill(rows)
}

# rbind data frames with differing columns (error rows lack the estimates).
.rbind_fill <- function(rows) {
  rows <- unlist(lapply(rows, function(r) if (is.data.frame(r)) split(r, seq_len(nrow(r))) else list()),
                 recursive = FALSE)
  cols <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(r) {
    r[setdiff(cols, names(r))] <- NA
    r[cols]
  }))
  rownames(out) <- NULL
  out
}
