\name{sim_mma_data}
\alias{sim_mma_data}
\alias{calibrate_eta0}
\title{Simulate mediation data with a partially missing mediator}
\description{
Generates \eqn{M = \alpha_A A + \alpha_X X + \epsilon_M},
\eqn{S_0 = \kappa_A A + \kappa_X X + \epsilon_S},
\eqn{Y = \beta_A A + \beta_M M + \beta_X X + \epsilon_Y} with
\eqn{corr(\epsilon_M, \epsilon_S) = \rho}; \eqn{S = S_0} or
\eqn{\exp(0.6 S_0)}. The mediator is MCAR, MAR given \eqn{(A, X)}, or MAR
given \eqn{(A, X, S)}; \code{calibrate_eta0} sets the missingness intercept to
reach a target missing rate. With \code{rho_Y}, an outcome surrogate
\eqn{S_Y = \kappa_X X + \epsilon_{SY}}, \eqn{corr(\epsilon_{SY}, \epsilon_Y) =}
\code{rho_Y}, is added and the outcome is made missing too; under MAR the
mediator's missingness depends on \eqn{(A, X, S, S_Y)} and the outcome's on
\eqn{(A, X, S, S_Y, M)}, with intercepts calibrated internally.
}
\usage{
sim_mma_data(n = 2000, alpha_A = 0.1, beta_M = 0.1, beta_A = 0.1,
             alpha_X = 0.3, beta_X = 0.3, kappa_A = 0, kappa_X = 0.3,
             rho = 0.6, miss_rate = 0.5, mechanism = c("MCAR", "MAR_AX", "MAR"),
             eta = c(A = 0.5, X = 0.5, S = 0.5),
             surrogate = c("linear", "nonlinear"),
             exposure = c("normal", "binary"), eta0 = NULL,
             rho_Y = NULL, miss_rate_Y = 0, mechanism_Y = c("MCAR", "MAR"),
             eta0_Y = NULL, errors = c("normal", "skewed"))
calibrate_eta0(miss_rate, mechanism = c("MAR", "MAR_AX"),
               eta = c(A = 0.5, X = 0.5, S = 0.5), n_ref = 2e5, ...)
}
\arguments{
\item{n}{sample size.}
\item{alpha_A, beta_M, beta_A, alpha_X, beta_X, kappa_A, kappa_X}{coefficients.}
\item{rho}{residual correlation between mediator and surrogate.}
\item{miss_rate}{marginal missing rate of the mediator.}
\item{mechanism}{missingness mechanism.}
\item{eta}{log odds ratios of observation in \eqn{A}, \eqn{X}, \eqn{S}.}
\item{surrogate}{\code{"linear"} or misspecified \code{"nonlinear"} surrogate.}
\item{exposure}{\code{"normal"} or \code{"binary"}.}
\item{eta0}{missingness intercept; calibrated when \code{NULL}.}
\item{rho_Y}{residual correlation between the outcome and its surrogate;
\code{NULL} for a complete outcome.}
\item{miss_rate_Y, mechanism_Y}{marginal missing rate and mechanism of the outcome.}
\item{eta0_Y}{intercept of the outcome missingness model; calibrated when \code{NULL}.}
\item{errors}{\code{"normal"} or \code{"skewed"} (standardized \eqn{\chi^2_4}) residuals.}
\item{n_ref}{reference sample size for calibration.}
\item{...}{data-generating parameters with the names used by \code{sim_mma_data}.}
}
\value{\code{sim_mma_data}: a data frame with \code{Y, A, X, S, M, M_full, R};
with \code{rho_Y}, also \code{S_Y, Y_full, R_Y}.
\code{calibrate_eta0}: the intercept.}
