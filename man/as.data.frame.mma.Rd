\name{as.data.frame.mma}
\alias{as.data.frame.mma}
\title{One-row summary of a mediation fit}
\description{
Returns a one-row data frame with the sample sizes, path estimates, effects,
tests (delta, joint significance and Monte Carlo p-values) and, for SynMediation, the decomposition of \eqn{\hat\beta_M} into its direct part
\eqn{\hat\omega_M} (\code{beta_direct}) and surrogate part
\eqn{\hat\nu^\top\hat\omega_S} (\code{beta_surrogate}), and the partial
\eqn{R^2} of the mediator on the surrogates given \eqn{Z}. For a binary
outcome fitted by sequential AIPW or IPW, also the missing-data ordering
(\code{missing_order}), the largest weight (\code{max_weight}), the effective
sample size (\code{ess}) and the largest single-subject share of the variance
of \eqn{\hat\beta_M} (\code{max_if_share_beta}).
}
\usage{
\method{as.data.frame}{mma}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
\item{x}{an object of class \code{"mma"}.}
\item{row.names, optional, ...}{ignored.}
}
\value{A one-row data frame.}
