\name{MissingMediatorAnalysis-package}
\alias{MissingMediatorAnalysis-package}
\alias{MissingMediatorAnalysis}
\docType{package}
\title{Mediation Analysis with a Partially Missing Mediator}
\description{
Product-of-coefficients mediation analysis when the mediator is missing at
random and a synthetic surrogate of the mediator is observed for everyone.
The exposure-mediator path is estimated by SynSurr (\code{\link{fit_synsurr}}),
the mediator-outcome path by the minimum-variance combination of complete-case
regression and the efficient AIPW estimator of Qin, Zhang and Leung (2017)
(\code{\link{fit_outcome_model}}), and \code{\link{mma}} combines
the two with stacked influence functions. For a binary outcome
(\code{outcome_family = "binomial"}), both paths are estimated by a
target-preserving sequential augmented inverse probability weighted estimator
whose validity depends only on the missingness models (see \code{\link{mma}}).
}
\references{
McCaw ZR, Gao J, Lin X, Gronsbell J (2024). Synthetic surrogates improve power
for genome-wide association studies of partially missing phenotypes in
population biobanks. \emph{Nature Genetics} 56, 1527--1536.

Qin J, Zhang B, Leung DHY (2017). Efficient augmented inverse probability
weighted estimation in missing data problems. \emph{Journal of Business &
Economic Statistics} 35, 86--97.
}
