\name{fit_synsurr}
\alias{fit_synsurr}
\title{Exposure-mediator path by synthetic surrogate (SynSurr) analysis}
\description{
Fits the bivariate model \eqn{M = Z^\top\alpha + \epsilon_M},
\eqn{S = Z^\top\alpha_S + \epsilon_S} with \eqn{M} partially missing. The
maximum likelihood estimate factors into two least-squares fits: \eqn{S} on
\eqn{Z} in all subjects, and \eqn{M} on \eqn{(S, Z)} in subjects with observed
\eqn{M}, giving \eqn{\hat\alpha = \hat\gamma + \hat\delta\hat\alpha_S}.
Standard errors are sandwich standard errors of the stacked estimating
equations. With \code{use}, the fit is restricted to a subset (rows outside
it have zero influence functions on the scale of all subjects); this is how the
outcome path with a missing outcome uses SynSurr among subjects with an
observed mediator.
}
\usage{fit_synsurr(M, S, Z, mediator_method = c("SynSurr", "CC"), use = NULL)}
\arguments{
\item{M}{mediator, \code{NA} where missing.}
\item{S}{synthetic surrogate, fully observed.}
\item{Z}{design matrix (with intercept) of exposure and covariates.}
\item{mediator_method}{\code{"SynSurr"} or complete-case \code{"CC"}.}
\item{use}{optional logical vector selecting the analysis population.}
}
\value{A list with \code{coef}, \code{vcov}, influence functions \code{IF}
(\eqn{n \times p}) and, for SynSurr, \code{alpha_S}, \code{delta}, \code{gamma}.}
