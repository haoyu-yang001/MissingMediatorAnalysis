#' Mediation analysis with a partially missing mediator and outcome
#'
#' Product-of-coefficients mediation analysis under the linear models
#' \deqn{M = \alpha_0 + \alpha_A A + X^\top\alpha_X + \epsilon_M,}
#' \deqn{Y = \beta_0 + \beta_A A + \beta_M M + X^\top\beta_X + \epsilon_Y,}
#' where the mediator \eqn{M}, and optionally the outcome \eqn{Y}, are
#' partially missing, a synthetic surrogate \eqn{S} of \eqn{M} is observed
#' for everyone, and, when \eqn{Y} is partially missing, so is a synthetic
#' surrogate of \eqn{Y} (\code{outcome_surrogate}).
#'
#' By default (\code{outcome_method = "SeqSynSurr"}) both paths come from
#' SynMediation: three least-squares fits, (1) each surrogate on
#' \eqn{(1, A, X)} in everyone, (2) \eqn{M} on the surrogates and
#' \eqn{(1, A, X)} among subjects with \eqn{M}, and (3) \eqn{Y} on
#' \eqn{M}, the surrogates and \eqn{(1, A, X)} among subjects with both
#' \eqn{M} and \eqn{Y}. \eqn{\alpha} (SynSurr) and \eqn{(\beta_A, \beta_M)}
#' follow in closed form by averaging over the surrogates. This is the maximum
#' likelihood estimator of the Gaussian model, equals ordinary least squares
#' with complete data, and requires that the mediator is missing at random
#' given \eqn{(A, X)} and the surrogates, and the outcome, among subjects with
#' \eqn{M}, given these and \eqn{M}. Subjects with \eqn{Y} but no \eqn{M}
#' are treated as missing \eqn{Y}.
#'
#' The natural indirect effect is \eqn{NIE = \alpha_A\beta_M}, the natural
#' direct effect \eqn{NDE = \beta_A}, and the total effect
#' \eqn{TE = NDE + NIE}. Both path estimators are asymptotically linear in the
#' same \eqn{n} subjects, so their influence functions are stacked and the
#' delta method gives a standard error for the NIE that includes
#' \eqn{Cov(\hat\alpha_A, \hat\beta_M)}. Two tests of \eqn{H_0: NIE = 0} are
#' reported: the delta-method (Sobel-type) Wald test and the joint-significance
#' test \eqn{p = \max(p_\alpha, p_\beta)}.
#'
#' Besides the delta-method (Wald) interval, a Monte Carlo confidence interval
#' for the NIE is reported. It draws \eqn{(\alpha_A^*, \beta_M^*)} from the
#' bivariate normal approximation to \eqn{(\hat\alpha_A, \hat\beta_M)},
#' including their estimated covariance, and takes quantiles of
#' \eqn{\alpha_A^*\beta_M^*}, so the interval reflects the skewness of a
#' product of two estimates. It uses the random number generator; call
#' \code{set.seed()} beforehand for reproducible limits.
#'
#' Alternatives for the mediator-outcome path (exposure-mediator path by
#' SynSurr): with a complete outcome, \code{"Combined"} (minimum-variance
#' combination of AIPW_new and complete case, valid when missingness does not
#' depend on \eqn{Y} given \eqn{(A, M, X)}) and \code{"AIPWnew"} (allows
#' dependence on \eqn{Y}); with a missing outcome, the doubly robust
#' \code{"SurrAIPW"}.
#'
#' @param data a data frame.
#' @param exposure,mediator,outcome,surrogate column names. The mediator may
#'   contain \code{NA}, and so may the outcome if \code{outcome_surrogate} is
#'   given; the other columns must be complete. For \code{"SeqSynSurr"},
#'   \code{surrogate} may name several columns: the synthetic surrogate of the
#'   mediator first, then auxiliary variables (fully observed, not in
#'   \code{covariates}) on which missingness may depend; they are integrated
#'   out like the surrogates.
#' @param outcome_surrogate optional column name of a synthetic surrogate of the
#'   outcome, observed for everyone. With a missing outcome it is required by
#'   \code{"SurrAIPW"} and \code{"SynSurr"}; for \code{"SeqSynSurr"} it is
#'   optional and weakens the MAR assumption for the outcome (missingness may
#'   then depend on it).
#' @param covariates character vector of covariate column names.
#' @param ps_vars variables (besides the intercept) for the propensity model
#'   and for the working model of \eqn{E[M \mid W]}. Defaults to
#'   \code{c(outcome, exposure, covariates, surrogate)}, or to
#'   \code{c(exposure, covariates, surrogate, outcome_surrogate)} when the
#'   outcome has missing values.
#' @param mediator_method \code{"SynSurr"} or \code{"CC"}.
#' @param outcome_method \code{"SeqSynSurr"} (default, see Details),
#'   \code{"Combined"}, \code{"SurrAIPW"}, \code{"SynSurr"}, \code{"AIPWnew"},
#'   \code{"AIPW"}, \code{"IPW"} or \code{"CC"} (see
#'   \code{\link{fit_outcome_model}}). With a missing outcome only
#'   \code{"SeqSynSurr"}, \code{"SurrAIPW"}, \code{"SynSurr"} and \code{"CC"}
#'   are available. \code{"SeqSynSurr"} requires
#'   \code{mediator_method = "SynSurr"}.
#' @param conf_level confidence level for intervals.
#' @param ps_trim optional lower bound for the fitted propensity.
#' @param mc_draws number of Monte Carlo draws for the NIE interval; 0 skips it.
#' @param na_action \code{"omit"} (default) drops, with a message, rows with
#'   missing values in the exposure, covariates, surrogates or \code{ps_vars};
#'   \code{"fail"} stops instead.
#' @param keep_IF keep the influence-function matrices in \code{fits}
#'   (\eqn{n \times p} each); set to \code{FALSE} to save memory when fitting
#'   many pathways.
#' @param verbose print messages about dropped rows.
#' @param basis for \code{"SeqSynSurr"}: \code{"linear"} (default) or a
#'   nonlinear model for \eqn{E(M \mid S, A, X)} in step 2, \code{"poly2"}
#'   (adds the square of each surrogate) or \code{"ns3"} (natural spline with
#'   3 df per surrogate); step 3 stays linear. The paths are then obtained by
#'   averaging the fitted moments over all subjects. Use it as a sensitivity
#'   analysis when \code{\link{mma_diagnostics}} flags nonlinearity in
#'   step 2. Only the average residual variance of step 2 enters the
#'   estimator, estimated among subjects with \eqn{M}; if this variance
#'   depends on the surrogates, check the \code{heteroskedasticity} part of
#'   \code{\link{mma_diagnostics}}. With a strongly heteroskedastic step 2 the
#'   nonlinear bases are not necessarily more reliable than the linear fit.
#' @param basis_outcome with \code{basis = "poly2"} or \code{"ns3"}, also use
#'   the basis of the surrogates in step 3 (sensitivity analysis; less
#'   efficient).
#' @return An object of class \code{"mma"}. Besides the estimates, the tests of
#'   \eqn{H_0: NIE = 0} include the p-value matching the Monte Carlo interval,
#'   \eqn{2\min\{P(\alpha^*\beta^* < 0), P(\alpha^*\beta^* > 0)\}} under the
#'   bivariate normal approximation, computed by numerical integration. Use
#'   \code{as.data.frame()} for a one-row summary.
#' @examples
#' dat <- sim_mma_data(n = 2000, miss_rate = 0.5, rho = 0.8)
#' fit <- mma(dat, exposure = "A", mediator = "M", outcome = "Y",
#'            surrogate = "S", covariates = "X")
#' fit
#' # Mediator and outcome both partially missing
#' dat2 <- sim_mma_data(n = 2000, miss_rate = 0.75, rho = 0.8,
#'                      rho_Y = 0.8, miss_rate_Y = 0.5)
#' mma(dat2, exposure = "A", mediator = "M", outcome = "Y", surrogate = "S",
#'     covariates = "X", outcome_surrogate = "S_Y")
#' @export
mma <- function(data, exposure, mediator, outcome, surrogate,
                covariates = NULL, ps_vars = NULL,
                mediator_method = c("SynSurr", "CC"),
                outcome_method = c("SeqSynSurr", "Combined", "SurrAIPW", "SynSurr",
                                   "AIPWnew", "AIPW", "IPW", "CC"),
                conf_level = 0.95, ps_trim = NULL, mc_draws = 1e5,
                outcome_surrogate = NULL, na_action = c("omit", "fail"),
                keep_IF = TRUE, verbose = TRUE, basis = c("linear", "poly2", "ns3"),
                basis_outcome = FALSE) {
  mediator_method <- match.arg(mediator_method)
  basis <- match.arg(basis)
  if (basis != "linear" && match.arg(outcome_method) != "SeqSynSurr")
    stop("basis is only used by outcome_method = \"SeqSynSurr\"")
  if (basis_outcome && basis == "linear") stop("basis_outcome = TRUE needs basis = \"poly2\" or \"ns3\"")
  data <- as.data.frame(data)
  y_missing <- anyNA(data[[outcome]])
  outcome_method <- match.arg(outcome_method)
  if (y_missing && !outcome_method %in% c("SeqSynSurr", "SurrAIPW", "SynSurr", "CC"))
    stop("with a missing outcome, outcome_method must be \"SeqSynSurr\", \"SurrAIPW\", \"SynSurr\" or \"CC\"")
  if (outcome_method == "SeqSynSurr" && mediator_method != "SynSurr")
    stop("outcome_method = \"SeqSynSurr\" requires mediator_method = \"SynSurr\"")
  if ((y_missing && !outcome_method %in% c("CC", "SeqSynSurr") || outcome_method == "SynSurr") &&
      is.null(outcome_surrogate))
    stop("with a missing outcome, outcome_method = \"SurrAIPW\" and \"SynSurr\" need outcome_surrogate")
  if (is.null(ps_vars)) ps_vars <- if (y_missing)
    c(exposure, covariates, surrogate, outcome_surrogate) else
    c(outcome, exposure, covariates, surrogate)
  if (y_missing && outcome %in% ps_vars) stop("ps_vars cannot contain a missing outcome")
  if (!is.numeric(data[[exposure]])) stop("the exposure must be numeric")
  complete_vars <- unique(c(exposure, surrogate, covariates, outcome_surrogate,
                            ps_vars, if (!y_missing) outcome))
  data <- .drop_incomplete(data, complete_vars, match.arg(na_action), verbose)

  n <- nrow(data)
  M <- data[[mediator]]
  Y <- data[[outcome]]
  Z <- .design(data, c(exposure, covariates))
  W <- .design(data, unique(ps_vars))
  S_Y <- if (is.null(outcome_surrogate)) NULL else data[[outcome_surrogate]]
  if (length(surrogate) > 1 && outcome_method != "SeqSynSurr")
    stop("several surrogate columns are only supported by outcome_method = \"SeqSynSurr\"")
  k_S <- length(surrogate) + !is.null(S_Y)
  .check_sizes(M, Y, ncol(Z), k_S * switch(basis, linear = 1, poly2 = 2, ns3 = 3), outcome_method)

  if (outcome_method == "SeqSynSurr") {
    # Both paths from one factorized fit that uses both surrogates.
    Zt <- cbind(Z, ifelse(is.na(M), 0, M))
    colnames(Zt) <- c(colnames(Z), mediator)
    S <- cbind(as.matrix(data[, surrogate, drop = FALSE]), S_Y)
    colnames(S) <- c(surrogate, outcome_surrogate)
    out <- .fit_seq_synsurr(Y, M, S, Z, Zt, basis, basis_outcome)
    med <- out$mediator
  } else {
    med <- fit_synsurr(M, data[[surrogate]], Z, mediator_method = mediator_method)
    out <- fit_outcome_model(Y, M, Z, W,
                             outcome_method = outcome_method,
                             mediator_name = mediator, ps_trim = ps_trim, S_Y = S_Y)
  }

  a <- med$coef[[exposure]]
  b <- out$coef[[mediator]]
  c_dir <- out$coef[[exposure]]
  IF_a <- med$IF[, exposure]
  IF_b <- out$IF[, mediator]
  IF_c <- out$IF[, exposure]

  nie <- a * b
  te <- c_dir + nie
  IF_nie <- b * IF_a + a * IF_b
  IF_te <- IF_c + IF_nie
  IF_pm <- (IF_nie * te - nie * IF_te) / te^2
  IF_eff <- cbind(alpha = IF_a, beta = IF_b, NIE = IF_nie, NDE = IF_c,
                  TE = IF_te, PM = IF_pm)
  V_eff <- .vcov_if(IF_eff)
  est <- c(alpha = a, beta = b, NIE = nie, NDE = c_dir, TE = te, PM = nie / te)
  effects <- .coef_table(est, sqrt(diag(V_eff)), conf_level)
  effects["PM", c("z", "p")] <- NA

  p_joint <- max(effects["alpha", "p"], effects["beta", "p"])
  p_mc <- .product_p(c(a, b), V_eff[c("alpha", "beta"), c("alpha", "beta")])
  tests <- data.frame(
    test = c("Delta (Sobel)", "Joint significance", "Monte Carlo"),
    statistic = c(effects["NIE", "z"], NA, NA),
    p = c(effects["NIE", "p"], p_joint, p_mc))

  nie_ci <- data.frame(method = "Delta", lower = effects["NIE", "lower"],
                       upper = effects["NIE", "upper"])
  if (mc_draws > 0) {
    mc <- .mc_product_ci(c(a, b), V_eff[c("alpha", "beta"), c("alpha", "beta")],
                         mc_draws, conf_level)
    nie_ci <- rbind(nie_ci, data.frame(method = "Monte Carlo", lower = mc[1],
                                       upper = mc[2]))
  }

  # Covariate coefficients that are not identified within an observed subset.
  med_tab <- .coef_table(med$coef, sqrt(diag(med$vcov)), conf_level)
  out_tab <- .coef_table(out$coef, sqrt(diag(out$vcov)), conf_level)
  med_tab[rownames(med_tab) %in% out$aliased, ] <- NA
  out_tab[rownames(out_tab) %in% out$aliased, ] <- NA
  if (!keep_IF) {
    med$IF <- NULL; out$IF <- NULL; out$components <- NULL
  }
  structure(list(
    effects = effects, tests = tests, nie_ci = nie_ci, conf_level = conf_level,
    vcov = V_eff,
    cor_ab = V_eff["alpha", "beta"] / sqrt(V_eff["alpha", "alpha"] * V_eff["beta", "beta"]),
    mediator_model = med_tab, outcome_model = out_tab,
    fits = list(mediator = med, outcome = out),
    methods = c(mediator = med$method, outcome = out$method),
    n = n, n_obs = sum(!is.na(M)), n_obs_Y = sum(!is.na(Y)),
    n_obs_both = sum(!is.na(M) & !is.na(Y)),
    n_Y_without_M = sum(is.na(M) & !is.na(Y)),
    aliased = out$aliased,
    vars = list(exposure = exposure, mediator = mediator, outcome = outcome,
                surrogate = surrogate, outcome_surrogate = outcome_surrogate,
                covariates = covariates, ps_vars = ps_vars),
    call = match.call()), class = "mma")
}

#' @export
print.mma <- function(x, digits = 4, ...) {
  cat("Mediation analysis with a partially missing mediator and outcome\n")
  cat(sprintf("  n = %d, mediator observed = %d (%.1f%% missing)\n",
              x$n, x$n_obs, 100 * (1 - x$n_obs / x$n)))
  if (!is.null(x$n_obs_Y) && x$n_obs_Y < x$n)
    cat(sprintf("  outcome observed = %d; mediator and outcome both observed = %d\n",
                x$n_obs_Y, x$n_obs_both))
  if (isTRUE(x$n_Y_without_M > 0) && x$n_obs_Y < x$n)
    cat(sprintf("  %d subjects with the outcome but no mediator were treated as missing the outcome\n",
                x$n_Y_without_M))
  if (length(x$aliased))
    cat("  not identified within the observed subset (coefficient set to NA):",
        paste(x$aliased, collapse = ", "), "\n")
  cat(sprintf("  exposure-mediator: %s;  mediator-outcome: %s\n\n",
              x$methods[["mediator"]], x$methods[["outcome"]]))
  print(signif(x$effects, digits))
  cat(sprintf("\n%g%% confidence intervals for the NIE\n", 100 * x$conf_level))
  print(x$nie_ci, digits = digits, row.names = FALSE)
  cat("\nTests of H0: NIE = 0\n")
  print(x$tests, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
summary.mma <- function(object, ...) {
  print(object, ...)
  cat("\nExposure-mediator model\n")
  print(signif(object$mediator_model, 4))
  cat("\nMediator-outcome model\n")
  print(signif(object$outcome_model, 4))
  cat(sprintf("\nEstimated corr(alpha_hat, beta_hat) = %.3f\n", object$cor_ab))
  invisible(object)
}

#' @export
coef.mma <- function(object, ...) {
  setNames(object$effects$estimate, rownames(object$effects))
}

#' @export
vcov.mma <- function(object, ...) object$vcov
