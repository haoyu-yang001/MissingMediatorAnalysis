\name{mma_diagnostics}
\alias{mma_diagnostics}
\alias{print.mma_diagnostics}
\title{Diagnostics for sequential SynSurr}
\description{
Checks the working assumptions of the default estimator for one pathway:
surrogate quality (partial \eqn{R^2} of the mediator on the surrogates given
\eqn{Z}; \eqn{R^2} of mediator and outcome on their surrogates); linearity
(robust Wald tests of squared terms added to steps 2 and 3); constant variance
in step 2 (Breusch-Pagan); the missingness mechanism (logistic models for
observing the mediator given \eqn{(Z, S)} and the outcome given
\eqn{(Z, S, M)}: Wald z, AUC and range of fitted probabilities); and a
Hausman-type comparison with complete-case analysis for \eqn{\alpha_A} and
\eqn{\beta_M}. A difference from complete case is expected when missingness
depends on the surrogates; without such dependence a large difference points to
misspecification. It also reports the change in \eqn{\hat\beta_M} implied
by a step-2 residual variance that depends on \eqn{(S, Z)} (the estimator only
uses its average among subjects with \eqn{M}; the average over all subjects is
estimated by inverse probability weighting; the change is about
\eqn{(\bar v_{all} - \bar v_{obs})(\omega_M - \beta_M) / Var(M \mid Z)}, assuming
a correctly specified step-2 mean), and
\eqn{\alpha_A}, \eqn{\beta_M} and the NIE refitted with a quadratic step 2
(\code{basis = "poly2"} in \code{\link{mma}}). All checks refer to the linear
fit, also when \code{fit} used a nonlinear basis.
}
\usage{
mma_diagnostics(fit, data)
\method{print}{mma_diagnostics}(x, digits = 3, ...)
}
\arguments{
\item{fit}{an \code{"mma"} object fitted with \code{outcome_method = "SeqSynSurr"}.}
\item{data}{the data frame used for \code{fit}.}
\item{x}{an object of class \code{"mma_diagnostics"}.}
\item{digits}{digits to print.}
\item{...}{ignored.}
}
\value{A list of data frames: \code{surrogate}, \code{linearity},
\code{variance}, \code{missingness}, \code{hausman},
\code{heteroskedasticity}, \code{nonlinearity}.}
\examples{
dat <- sim_mma_data(n = 2000, miss_rate = 0.75, rho = 0.8,
                    rho_Y = 0.8, miss_rate_Y = 0.5)
fit <- mma(dat, "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y")
mma_diagnostics(fit, dat)
}
