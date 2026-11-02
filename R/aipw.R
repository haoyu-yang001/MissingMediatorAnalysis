#' Mediator-outcome path with a partially missing mediator
#'
#' Estimates the outcome regression
#' \deqn{Y = C^\top\beta_C + \beta_M M + \epsilon_Y}
#' when the regressor \eqn{M} is missing at random given the fully observed
#' vector \eqn{W} (typically \eqn{W = (1, Y, A, X, S)}). The full-data
#' estimating function is \eqn{U(\beta) = \tilde Z (Y - \tilde Z^\top\beta)}
#' with \eqn{\tilde Z = (C, M)}.
#'
#' \describe{
#'   \item{\code{"CC"}}{complete-case least squares.}
#'   \item{\code{"IPW"}}{inverse probability weighting,
#'     \eqn{E_n[R\,U(\beta)/\omega(W)] = 0}.}
#'   \item{\code{"AIPW"}}{Robins, Rotnitzky and Zhao (1994):
#'     \eqn{E_n[R\,U/\omega - (R-\omega)\,q(W;\beta)/\omega] = 0} with
#'     \eqn{q(W;\beta)} a working model for \eqn{E[U(\beta)\mid W]}.}
#'   \item{\code{"AIPWnew"}}{Qin, Zhang and Leung (2017): additionally subtracts
#'     \eqn{\hat\tau A(W; \hat q, \tilde q)}, the least-squares projection of the
#'     AIPW estimating function on the span of the augmentation term and the
#'     propensity score, with \eqn{\hat\tau = \hat B_1 \hat B_2^{-1}}. Iterated
#'     to convergence.}
#'   \item{\code{"Combined"}}{minimum-variance linear combination of
#'     \code{"AIPWnew"} and \code{"CC"},
#'     \eqn{\hat\beta_{new} - \Gamma(\hat\beta_{new} - \hat\beta_{CC})}, with
#'     \eqn{\Gamma} estimated from the influence functions. CC is consistent
#'     only if missingness does not depend on \eqn{Y} given \eqn{(C, M)}, so
#'     this method requires that assumption; it is then asymptotically at
#'     least as efficient as either component.}
#'   \item{\code{"SurrAIPW"}}{for a partially missing outcome (default in
#'     \code{\link{mma}} in that case): surrogate-augmented AIPW for monotone
#'     missingness. With \eqn{R_1} the indicator that \eqn{M} is observed and
#'     \eqn{R_2} that both \eqn{M} and \eqn{Y} are, it solves
#'     \eqn{E_n[R_2 U/\pi_2 + (1 - R_1/\pi_1) q_1(W) + (R_1/\pi_1 - R_2/\pi_2)
#'     q_2(W, M)] = 0}, where \eqn{\pi_1 = P(R_1 = 1 \mid W)},
#'     \eqn{\pi_2 = \pi_1 P(R_2 = 1 \mid R_1 = 1, W, M)} (both logistic),
#'     \eqn{q_1 = E[U \mid W]} and \eqn{q_2 = E[U \mid W, M]}. \eqn{W} should
#'     contain the surrogates of both \eqn{M} and \eqn{Y}; \eqn{q_2} uses a
#'     linear working model for \eqn{Y} given \eqn{(W, M)}, the SynSurr idea of
#'     predicting the missing outcome from its surrogate given the regressors.
#'     Requires that \eqn{M} is missing at random given \eqn{W} and \eqn{Y}
#'     given \eqn{(W, M)}. Subjects with \eqn{Y} but no \eqn{M} are treated as
#'     missing both.}
#'   \item{\code{"SynSurr"}}{for a partially missing outcome: SynSurr
#'     (\code{\link{fit_synsurr}}) with target \eqn{Y}, surrogate \code{S_Y}
#'     and regressors \eqn{\tilde Z}, among subjects with an observed mediator.
#'     Requires that \eqn{Y} is missing at random given
#'     \eqn{(\tilde Z, S_Y)} and that missingness of \eqn{M} does not depend on
#'     \eqn{Y} given \eqn{\tilde Z}.}
#' }
#'
#' Only \code{"SurrAIPW"}, \code{"SynSurr"} and \code{"CC"} allow a missing
#' outcome (see also \code{outcome_method = "SeqSynSurr"} in \code{\link{mma}});
#' \code{"CC"}
#' then uses subjects with both \eqn{M} and \eqn{Y} observed.
#'
#' The propensity \eqn{\omega(W) = P(R = 1 \mid W)} is logistic in \eqn{W}. The
#' working model for \eqn{q} is \eqn{M \mid W \sim N(W^\top\gamma, \sigma^2)},
#' fitted on complete cases, which gives \eqn{E[M\mid W]} and
#' \eqn{E[M^2\mid W]} and hence \eqn{E[U(\beta)\mid W]} in closed form.
#' Standard errors are sandwich standard errors that account for estimating
#' the propensity score.
#'
#' @param Y numeric outcome; may contain \code{NA} only for
#'   \code{"SurrAIPW"}, \code{"SynSurr"} and \code{"CC"}.
#' @param M numeric mediator, \code{NA} where missing.
#' @param C design matrix (with intercept) of the fully observed regressors of
#'   the outcome model (exposure and covariates).
#' @param W design matrix (with intercept) of the fully observed variables used
#'   by the propensity and working models.
#' @param outcome_method one of \code{"Combined"} (default), \code{"SurrAIPW"},
#'   \code{"SynSurr"},
#'   \code{"AIPWnew"}, \code{"AIPW"}, \code{"IPW"}, \code{"CC"}.
#' @param S_Y synthetic surrogate of the outcome, required for
#'   \code{"SynSurr"}; must be observed for subjects with an observed mediator.
#' @param mediator_name column name for the mediator coefficient.
#' @param ps_trim optional lower bound for the fitted propensity.
#' @param maxit,tol iteration control for \code{"AIPWnew"}.
#' @return A list with \code{coef}, \code{vcov}, the influence-function matrix
#'   \code{IF}, \code{method}, and for weighted methods the propensity
#'   coefficients \code{eta} and fitted \code{omega}. For \code{"Combined"},
#'   \code{weight_cc} (the matrix \eqn{\Gamma}) and the two component fits
#'   \code{components}.
#' @export
fit_outcome_model <- function(Y, M, C, W = NULL,
                              outcome_method = c("Combined", "SurrAIPW", "SynSurr",
                                                 "AIPWnew", "AIPW", "IPW", "CC"),
                              mediator_name = "M", ps_trim = NULL,
                              maxit = 100, tol = 1e-10, S_Y = NULL) {
  outcome_method <- match.arg(outcome_method)
  C <- as.matrix(C)
  if (is.null(colnames(C))) colnames(C) <- paste0("C", seq_len(ncol(C)))
  n <- length(Y)
  R <- as.numeric(!is.na(M))
  Zt <- cbind(C, ifelse(R == 1, M, 0))
  colnames(Zt) <- c(colnames(C), mediator_name)
  K <- ncol(Zt)
  pc <- ncol(C)

  if (outcome_method == "SynSurr") {
    if (is.null(S_Y)) stop("outcome_method = \"SynSurr\" requires S_Y")
    if (anyNA(S_Y[R == 1])) stop("S_Y must be observed when the mediator is")
    fit <- fit_synsurr(Y, S_Y, Zt, use = R == 1)
    fit$n_obs_Y <- sum(R == 1 & !is.na(Y))
    return(fit)
  }
  if (outcome_method == "SurrAIPW") {
    if (is.null(W) || anyNA(W)) stop("SurrAIPW requires a complete W")
    return(.fit_surr_aipw(Y, M, C, as.matrix(W), Zt, ps_trim))
  }
  if (anyNA(Y)) {
    if (outcome_method != "CC")
      stop("a missing outcome requires outcome_method = \"SurrAIPW\", \"SynSurr\" or \"CC\"")
    fit <- .wls_if(Zt, Y, w = R * !is.na(Y))
    return(list(coef = fit$coef, vcov = .vcov_if(fit$IF), IF = fit$IF,
                method = "CC"))
  }
  W <- as.matrix(W)

  if (outcome_method == "Combined" && any(R == 0)) {
    new <- fit_outcome_model(Y, M, C, W, "AIPWnew", mediator_name, ps_trim,
                             maxit, tol)
    cc <- fit_outcome_model(Y, M, C, W, "CC", mediator_name)
    comb <- .combine_if(new$coef, new$IF, cc$coef, cc$IF)
    return(list(coef = comb$coef, vcov = .vcov_if(comb$IF), IF = comb$IF,
                method = "Combined", weight_cc = comb$Gamma,
                components = list(AIPWnew = new, CC = cc)))
  }

  if (outcome_method %in% c("CC", "Combined") || all(R == 1)) {
    fit <- .wls_if(Zt, Y, w = R)
    return(list(coef = fit$coef, vcov = .vcov_if(fit$IF), IF = fit$IF,
                method = if (all(R == 1)) "full" else "CC"))
  }

  ps <- .fit_propensity(W, R, ps_trim)
  omega <- ps$omega
  wR <- R / omega
  a <- (R - omega) / omega

  if (outcome_method == "IPW") {
    fit <- .wls_if(Zt, Y, w = wR)
    psi <- .residualize(fit$IF %*% fit$J, ps$score)
    IF <- psi %*% .pinv(fit$J)
    colnames(IF) <- colnames(Zt)
    return(list(coef = fit$coef, vcov = .vcov_if(IF), IF = IF, method = "IPW",
                eta = ps$eta, omega = omega))
  }

  # Working model M | W ~ N(W'gamma, sigma^2) from complete cases.
  wm <- .wls_if(W, M, w = R)
  m1 <- drop(W %*% wm$coef)
  m2 <- m1^2 + sum(wm$resid[R == 1]^2) / (sum(R) - ncol(W))

  # q(W; beta) = E[Zt Y | W] - E[Zt Zt' | W] beta, linear in beta.
  EZY <- cbind(C * Y, m1 * Y)
  qfun <- function(beta) {
    lc <- drop(C %*% beta[1:pc])
    EZY - cbind(C * (lc + m1 * beta[K]), m1 * lc + m2 * beta[K])
  }
  # sum_i c_i E[Zt_i Zt_i' | W_i]
  sum_EZZ <- function(c) {
    cm <- crossprod(C, c * m1)
    rbind(cbind(crossprod(C * c, C), cm), c(cm, sum(c * m2)))
  }
  Ufun <- function(beta) Zt * (R * drop(Y - Zt %*% beta))

  # AIPW: sum_i [wR_i U_i(beta) - a_i q_i(beta)] = 0 is linear in beta.
  A_mat <- crossprod(Zt * wR, Zt) - sum_EZZ(a)
  b_vec <- drop(crossprod(Zt * wR, Y) - crossprod(EZY, a))
  beta <- drop(solve(A_mat, b_vec))
  n_iter <- 0L

  if (outcome_method == "AIPW") {
    psi <- wR * Ufun(beta) - a * qfun(beta)
  } else {
    # AIPW_new: q-tilde = omega_eta / (1 - omega) = omega * W for the logit link.
    q_tilde <- omega * W
    for (it in seq_len(maxit)) {
      U <- Ufun(beta)
      q <- qfun(beta)
      h <- cbind(q, q_tilde)
      B1 <- crossprod((U - q) * (R * (1 - omega) / omega^2), h) / n
      B2 <- crossprod(h * a) / n
      tau <- B1 %*% .pinv(B2)
      beta_new <- drop(solve(A_mat, b_vec - tau %*% colSums(h * a)))
      n_iter <- it
      done <- max(abs(beta_new - beta)) < tol * (1 + max(abs(beta)))
      beta <- beta_new
      if (done) break
    }
    if (n_iter == maxit) warning("AIPWnew did not converge")
    U <- Ufun(beta)
    q <- qfun(beta)
    h <- cbind(q, q_tilde)
    psi <- wR * U - a * q - (h * a) %*% t(tau)
  }

  names(beta) <- colnames(Zt)
  J <- A_mat / n
  IF <- .residualize(psi, ps$score) %*% t(solve(J))
  colnames(IF) <- colnames(Zt)
  list(coef = beta, vcov = .vcov_if(IF), IF = IF, method = outcome_method,
       eta = ps$eta, omega = omega, working_model = wm$coef, n_iter = n_iter)
}
