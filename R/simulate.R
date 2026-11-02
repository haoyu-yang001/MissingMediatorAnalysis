#' Simulate mediation data with a partially missing mediator
#'
#' Generates
#' \deqn{M = \alpha_A A + \alpha_X X + \epsilon_M,\quad
#'       S_0 = \kappa_A A + \kappa_X X + \epsilon_S,\quad
#'       Y = \beta_A A + \beta_M M + \beta_X X + \epsilon_Y,}
#' with \eqn{X, \epsilon_Y \sim N(0,1)} and \eqn{(\epsilon_M, \epsilon_S)}
#' standard bivariate normal with correlation \code{rho}. The synthetic
#' surrogate is \eqn{S = S_0} (\code{surrogate = "linear"}) or
#' \eqn{S = \exp(0.6 S_0)} (\code{"nonlinear"}, a misspecified surrogate).
#' Missingness of \eqn{M}:
#' \describe{
#'   \item{\code{"MCAR"}}{\eqn{P(R = 1) = 1 - } \code{miss_rate};}
#'   \item{\code{"MAR_AX"}}{\eqn{logit P(R=1) = \eta_0 + \eta_A A + \eta_X X};}
#'   \item{\code{"MAR"}}{\eqn{logit P(R=1) = \eta_0 + \eta_A A + \eta_X X + \eta_S S};}
#' }
#' with \eqn{\eta_0} calibrated so the marginal missing rate is
#' \code{miss_rate}.
#'
#' If \code{rho_Y} is given, an outcome surrogate
#' \eqn{S_Y = \kappa_X X + \epsilon_{SY}} with
#' \eqn{\mathrm{corr}(\epsilon_{SY}, \epsilon_Y) = } \code{rho_Y} is added
#' (an earlier measurement of the outcome, say) and the outcome is also made
#' missing, with marginal rate \code{miss_rate_Y}: completely at random, or
#' (\code{mechanism_Y = "MAR"}) with
#' \eqn{logit P(R_Y=1) = \eta_{0Y} + 0.5(A + X + S + S_Y + M)}, which is missing
#' at random given the fully observed variables and the mediator. With
#' \code{rho_Y} and \code{mechanism = "MAR"}, missingness of the mediator also
#' depends on the outcome surrogate:
#' \eqn{logit P(R=1) = \eta_0 + 0.5(A + X + S + S_Y)}, and both intercepts are
#' calibrated internally unless given.
#'
#' @param n sample size.
#' @param alpha_A,beta_M,beta_A path coefficients; the true NIE is
#'   \code{alpha_A * beta_M}.
#' @param alpha_X,beta_X,kappa_A,kappa_X remaining coefficients.
#' @param rho residual correlation between mediator and surrogate.
#' @param miss_rate marginal probability that \eqn{M} is missing.
#' @param mechanism \code{"MCAR"}, \code{"MAR_AX"} or \code{"MAR"}.
#' @param eta named vector \code{c(A =, X =, S =)} of missingness log odds
#'   ratios.
#' @param surrogate \code{"linear"} or \code{"nonlinear"}.
#' @param exposure \code{"normal"} or \code{"binary"} (Bernoulli(0.5)).
#' @param eta0 intercept of the missingness model; calibrated by
#'   \code{\link{calibrate_eta0}} when \code{NULL}.
#' @param rho_Y residual correlation between the outcome and its surrogate;
#'   \code{NULL} (default) for a fully observed outcome.
#' @param miss_rate_Y,mechanism_Y marginal missing rate and mechanism
#'   (\code{"MCAR"} or \code{"MAR"}) of the outcome.
#' @param eta0_Y intercept of the outcome missingness model (calibrated when
#'   \code{NULL}).
#' @param errors \code{"normal"} or \code{"skewed"} (standardized
#'   \eqn{\chi^2_4}) residuals for \eqn{M}, \eqn{Y} and the surrogates.
#' @return A data frame with columns \code{Y, A, X, S, M} (with \code{NA}),
#'   \code{M_full} and \code{R}; with \code{rho_Y}, also \code{S_Y},
#'   \code{Y_full} and \code{R_Y}, and \code{Y} has \code{NA}.
#' @export
sim_mma_data <- function(n = 2000, alpha_A = 0.1, beta_M = 0.1, beta_A = 0.1,
                         alpha_X = 0.3, beta_X = 0.3, kappa_A = 0,
                         kappa_X = 0.3, rho = 0.6, miss_rate = 0.5,
                         mechanism = c("MCAR", "MAR_AX", "MAR"),
                         eta = c(A = 0.5, X = 0.5, S = 0.5),
                         surrogate = c("linear", "nonlinear"),
                         exposure = c("normal", "binary"), eta0 = NULL,
                         rho_Y = NULL, miss_rate_Y = 0,
                         mechanism_Y = c("MCAR", "MAR"), eta0_Y = NULL,
                         errors = c("normal", "skewed")) {
  mechanism <- match.arg(mechanism)
  surrogate <- match.arg(surrogate)
  exposure <- match.arg(exposure)
  mechanism_Y <- match.arg(mechanism_Y)
  errors <- match.arg(errors)
  gen <- function(n) .gen_full(n, alpha_A, beta_M, beta_A, alpha_X, beta_X, kappa_A,
                               kappa_X, rho, surrogate, exposure, rho_Y, errors)
  if (!is.null(rho_Y))
    return(.sim_ymiss(gen(n), gen, miss_rate, mechanism, miss_rate_Y, mechanism_Y,
                      eta0, eta0_Y))
  dat <- gen(n)
  lin <- .miss_linpred(dat, mechanism, eta)
  if (mechanism == "MCAR") {
    p_obs <- rep(1 - miss_rate, n)
  } else {
    if (is.null(eta0)) {
      eta0 <- calibrate_eta0(miss_rate, mechanism, eta, alpha_A = alpha_A,
                             beta_M = beta_M, beta_A = beta_A,
                             alpha_X = alpha_X, beta_X = beta_X,
                             kappa_A = kappa_A, kappa_X = kappa_X, rho = rho,
                             surrogate = surrogate, exposure = exposure)
    }
    p_obs <- plogis(eta0 + lin)
  }
  dat$R <- rbinom(n, 1, p_obs)
  dat$M <- ifelse(dat$R == 1, dat$M_full, NA)
  dat[, c("Y", "A", "X", "S", "M", "M_full", "R")]
}

# Missing mediator and outcome. Under MAR the mediator's missingness depends on
# (A, X, S, S_Y) and the outcome's on (A, X, S, S_Y, M).
.sim_ymiss <- function(dat, gen, miss_rate, mechanism, miss_rate_Y, mechanism_Y,
                       eta0 = NULL, eta0_Y = NULL) {
  n <- nrow(dat)
  lin_M <- function(d) 0.5 * (d$A + d$X + d$S + d$S_Y)
  lin_Y <- function(d) 0.5 * (d$A + d$X + d$S + d$S_Y + d$M_full)
  if ((mechanism != "MCAR" && is.null(eta0)) || (mechanism_Y == "MAR" && is.null(eta0_Y))) {
    e0 <- .calib_ymiss(gen, miss_rate, miss_rate_Y)
    if (is.null(eta0)) eta0 <- e0[[1]]
    if (is.null(eta0_Y)) eta0_Y <- e0[[2]]
  }
  p_M <- if (mechanism == "MCAR") rep(1 - miss_rate, n) else plogis(eta0 + lin_M(dat))
  p_Y <- if (mechanism_Y == "MCAR") rep(1 - miss_rate_Y, n) else plogis(eta0_Y + lin_Y(dat))
  dat$R <- rbinom(n, 1, p_M)
  dat$R_Y <- rbinom(n, 1, p_Y)
  dat$M <- ifelse(dat$R == 1, dat$M_full, NA)
  dat$Y_full <- dat$Y
  dat$Y <- ifelse(dat$R_Y == 1, dat$Y_full, NA)
  dat[, c("Y", "A", "X", "S", "S_Y", "M", "M_full", "Y_full", "R", "R_Y")]
}

# Intercepts of the MAR models in .sim_ymiss giving the target missing rates.
.calib_ymiss <- function(gen, miss_rate, miss_rate_Y, n_ref = 2e5) {
  ref <- gen(n_ref)
  root <- function(lin, rate) uniroot(function(e) mean(plogis(e + lin)) - (1 - rate),
                                      c(-30, 30), tol = 1e-10)$root
  c(root(0.5 * (ref$A + ref$X + ref$S + ref$S_Y), miss_rate),
    root(0.5 * (ref$A + ref$X + ref$S + ref$S_Y + ref$M_full), miss_rate_Y))
}

#' Calibrate the missingness intercept
#'
#' Finds \eqn{\eta_0} such that the marginal missing rate equals
#' \code{miss_rate}, using a large reference sample from the same
#' data-generating model.
#'
#' @inheritParams sim_mma_data
#' @param n_ref size of the reference sample.
#' @param ... data-generating parameters passed on to the internal generator
#'   (same names as in \code{\link{sim_mma_data}}).
#' @return The calibrated intercept.
#' @export
calibrate_eta0 <- function(miss_rate, mechanism = c("MAR", "MAR_AX"),
                           eta = c(A = 0.5, X = 0.5, S = 0.5),
                           n_ref = 2e5, ...) {
  mechanism <- match.arg(mechanism)
  args <- modifyList(
    list(n = n_ref, alpha_A = 0.1, beta_M = 0.1, beta_A = 0.1, alpha_X = 0.3,
         beta_X = 0.3, kappa_A = 0, kappa_X = 0.3, rho = 0.6,
         surrogate = "linear", exposure = "normal"), list(...))
  ref <- do.call(.gen_full, args)
  lin <- .miss_linpred(ref, mechanism, eta)
  uniroot(function(e0) mean(plogis(e0 + lin)) - (1 - miss_rate),
          c(-30, 30), tol = 1e-10)$root
}

.gen_full <- function(n, alpha_A, beta_M, beta_A, alpha_X, beta_X, kappa_A,
                      kappa_X, rho, surrogate, exposure, rho_Y = NULL,
                      errors = "normal") {
  rerr <- if (errors == "normal") rnorm else function(n) (rchisq(n, 4) - 4) / sqrt(8)
  A <- if (exposure == "normal") rnorm(n) else rbinom(n, 1, 0.5)
  X <- rnorm(n)
  e_M <- rerr(n)
  e_S <- rho * e_M + sqrt(1 - rho^2) * rerr(n)
  M <- alpha_A * A + alpha_X * X + e_M
  S <- kappa_A * A + kappa_X * X + e_S
  if (surrogate == "nonlinear") S <- exp(0.6 * S)
  e_Y <- rerr(n)
  Y <- beta_A * A + beta_M * M + beta_X * X + e_Y
  out <- data.frame(Y = Y, A = A, X = X, S = S, M_full = M)
  if (!is.null(rho_Y))
    out$S_Y <- kappa_X * X + rho_Y * e_Y + sqrt(1 - rho_Y^2) * rerr(n)
  out
}

.miss_linpred <- function(dat, mechanism, eta) {
  switch(mechanism,
         MCAR = rep(0, nrow(dat)),
         MAR_AX = eta[["A"]] * dat$A + eta[["X"]] * dat$X,
         MAR = eta[["A"]] * dat$A + eta[["X"]] * dat$X +
           eta[["S"]] * dat$S)
}
