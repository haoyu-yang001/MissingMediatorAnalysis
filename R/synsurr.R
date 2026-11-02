#' Exposure-mediator path by synthetic surrogate (SynSurr) analysis
#'
#' Fits the bivariate model
#' \deqn{M = Z^\top\alpha + \epsilon_M,\quad S = Z^\top\alpha_S + \epsilon_S,}
#' where the mediator \eqn{M} is partially missing and the synthetic surrogate
#' \eqn{S} is observed for everyone (McCaw et al., 2024, Nature Genetics).
#' Under bivariate normality the maximum likelihood estimate factors into two
#' least-squares fits: (1) \eqn{S} on \eqn{Z} in all \eqn{n} subjects, giving
#' \eqn{\hat\alpha_S}; (2) \eqn{M} on \eqn{(S, Z)} in the subjects with
#' observed \eqn{M}, giving \eqn{(\hat\delta, \hat\gamma)}. Then
#' \eqn{\hat\alpha = \hat\gamma + \hat\delta\hat\alpha_S}. Standard errors come
#' from the sandwich variance of the stacked estimating equations.
#'
#' The same function fits the mediator-outcome path when the outcome is
#' partially missing: the target is then \eqn{Y}, the surrogate a synthetic
#' surrogate of \eqn{Y}, \eqn{Z = (1, A, X, M)}, and \code{use} restricts the
#' fit to subjects with an observed mediator.
#'
#' @param M numeric mediator, \code{NA} where missing.
#' @param S numeric synthetic surrogate, fully observed.
#' @param Z numeric design matrix (with intercept) of the exposure and
#'   covariates, fully observed.
#' @param mediator_method \code{"SynSurr"} (default) or \code{"CC"} for the
#'   complete-case regression of \code{M} on \code{Z}.
#' @param use optional logical vector selecting the subjects that form the
#'   analysis population (default: all). Rows outside it do not contribute and
#'   may contain \code{NA}; their influence functions are zero, and the
#'   influence functions stay on the scale of all \code{length(M)} subjects so
#'   that they can be stacked with fits on the full sample.
#' @return A list with \code{coef}, \code{vcov}, the \eqn{n \times p}
#'   influence-function matrix \code{IF}, and for SynSurr the components
#'   \code{alpha_S}, \code{delta}, \code{gamma}.
#' @export
fit_synsurr <- function(M, S, Z, mediator_method = c("SynSurr", "CC"),
                        use = NULL) {
  mediator_method <- match.arg(mediator_method)
  Z <- as.matrix(Z)
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  u <- if (is.null(use)) rep(1, length(M)) else as.numeric(use)
  R <- as.numeric(!is.na(M)) * u

  if (mediator_method == "CC") {
    fit <- .wls_if(Z, M, w = R)
    return(list(coef = fit$coef, vcov = .vcov_if(fit$IF), IF = fit$IF,
                method = "CC"))
  }

  # Step 1: surrogate on Z in all subjects of the analysis population.
  fit_S <- .wls_if(Z, S, w = u)
  # Step 2: mediator on (S, Z) among subjects with observed mediator.
  X2 <- cbind(S = S, Z)
  fit_M <- .wls_if(X2, M, w = R)
  delta <- fit_M$coef[1]
  gamma <- fit_M$coef[-1]

  coef <- gamma + delta * fit_S$coef
  IF <- fit_M$IF[, -1, drop = FALSE] + delta * fit_S$IF +
    outer(fit_M$IF[, 1], fit_S$coef)
  names(coef) <- colnames(IF) <- colnames(Z)

  list(coef = coef, vcov = .vcov_if(IF), IF = IF, method = "SynSurr",
       alpha_S = fit_S$coef, delta = unname(delta), gamma = gamma)
}
