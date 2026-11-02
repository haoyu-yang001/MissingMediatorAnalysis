\name{fit_outcome_model}
\alias{fit_outcome_model}
\title{Mediator-outcome path with a partially missing mediator}
\description{
Estimates \eqn{Y = C^\top\beta_C + \beta_M M + \epsilon_Y} when \eqn{M} is
missing at random given the fully observed \eqn{W}. Methods: complete case,
IPW, AIPW (Robins, Rotnitzky and Zhao, 1994) and AIPW_new (Qin, Zhang and
Leung, 2017), which subtracts from the AIPW estimating function its
least-squares projection on the augmentation term and the propensity score.
\code{"Combined"} is the minimum-variance linear combination of AIPW_new and
complete case, \eqn{\hat\beta_{new} - \Gamma(\hat\beta_{new} - \hat\beta_{CC})},
with \eqn{\Gamma} estimated from the influence functions; it requires that
missingness does not depend on \eqn{Y} given \eqn{(C, M)}, under which CC is
consistent.
The propensity is logistic in \eqn{W}; the working model is
\eqn{M \mid W \sim N(W^\top\gamma, \sigma^2)}, which gives
\eqn{E[U(\beta) \mid W]} in closed form.

For a partially missing outcome: \code{"SurrAIPW"} solves the monotone-coarsening
AIPW equation \eqn{E_n[R_2 U/\pi_2 + (1 - R_1/\pi_1) q_1(W) + (R_1/\pi_1 -
R_2/\pi_2) q_2(W, M)] = 0}, where \eqn{R_1} indicates an observed mediator,
\eqn{R_2} an observed mediator and outcome, \eqn{\pi_1 = P(R_1 = 1 \mid W)},
\eqn{\pi_2 = \pi_1 P(R_2 = 1 \mid R_1 = 1, W, M)}, \eqn{q_1 = E[U \mid W]}
and \eqn{q_2 = E[U \mid W, M]}; \eqn{W} should contain the surrogates of both
\eqn{M} and \eqn{Y}. \code{"SynSurr"} applies \code{\link{fit_synsurr}} to
\eqn{Y} with surrogate \code{S_Y} among subjects with an observed mediator.
\code{"CC"} then uses subjects with both observed.
}
\usage{
fit_outcome_model(Y, M, C, W = NULL,
                  outcome_method = c("Combined", "SurrAIPW", "SynSurr", "AIPWnew",
                                     "AIPW", "IPW", "CC"),
                  mediator_name = "M", ps_trim = NULL, maxit = 100, tol = 1e-10,
                  S_Y = NULL)
}
\arguments{
\item{Y}{outcome; may contain \code{NA} only for \code{"SurrAIPW"}, \code{"SynSurr"} and \code{"CC"}.}
\item{M}{mediator, \code{NA} where missing.}
\item{C}{design matrix (with intercept) of the fully observed regressors.}
\item{W}{design matrix (with intercept) for the propensity and working models.}
\item{outcome_method}{estimator: \code{"Combined"} (default), \code{"SurrAIPW"}, \code{"SynSurr"}, \code{"AIPWnew"},
\code{"AIPW"}, \code{"IPW"} or \code{"CC"}.}
\item{mediator_name}{name of the mediator coefficient.}
\item{ps_trim}{optional lower bound for fitted propensities.}
\item{maxit, tol}{iteration control for AIPW_new.}
\item{S_Y}{synthetic surrogate of the outcome, required for \code{"SynSurr"}.}
}
\value{A list with \code{coef}, \code{vcov}, influence functions \code{IF},
\code{method}, and for weighted methods \code{eta} and \code{omega}.
For \code{"Combined"}, \code{weight_cc} (the matrix \eqn{\Gamma}) and the two
component fits \code{components}.}
