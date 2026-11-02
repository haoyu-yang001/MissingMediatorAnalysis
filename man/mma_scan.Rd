\name{mma_scan}
\alias{mma_scan}
\title{Mediation analysis over many exposures and mediators}
\description{
Calls \code{\link{mma}} for every exposure x mediator pair and returns one row
per pair. Errors and warnings are recorded in the columns \code{status} and
\code{message} rather than stopping the scan. Rows with a missing exposure are
dropped pair by pair.
}
\usage{
mma_scan(data, exposures, mediators, surrogates, outcome, covariates = NULL,
         outcome_surrogate = NULL, auxiliary = NULL, ..., n_cores = 1,
         seed = NULL)
}
\arguments{
\item{data}{a data frame with the mediators, outcome, surrogates and covariates.}
\item{exposures}{a character vector of exposure columns of \code{data}, or a
numeric matrix with \code{nrow(data)} rows and column names (e.g. thousands of
proteins).}
\item{mediators}{character vector of mediator columns.}
\item{surrogates}{character vector of surrogate columns, one per mediator.}
\item{outcome, covariates, outcome_surrogate}{as in \code{\link{mma}}.}
\item{auxiliary}{optional fully observed auxiliary columns added to every
mediator's surrogate (sequential SynSurr only; see \code{\link{mma}}).}
\item{...}{further arguments to \code{\link{mma}}, e.g. \code{outcome_method},
\code{mediator_method}, \code{mc_draws}.}
\item{n_cores}{number of forked cores (\code{parallel::mclapply}). Exposures
whose worker fails (e.g. out of memory) are returned as error rows.}
\item{seed}{optional integer; if given, the Monte Carlo interval for the
\eqn{j}-th exposure uses seed \code{seed + j}, so results are reproducible also
in parallel.}
}
\value{
A data frame with columns \code{exposure}, \code{mediator}, \code{status}
(\code{"ok"}, \code{"warning"} or \code{"error"}), \code{message} and the
columns of \code{\link{as.data.frame.mma}}.
}
\examples{
dat <- sim_mma_data(n = 2000, miss_rate = 0.75, rho = 0.8,
                    rho_Y = 0.8, miss_rate_Y = 0.5)
G <- cbind(A = dat$A, A2 = rnorm(2000))
mma_scan(dat, G, mediators = "M", surrogates = "S", outcome = "Y",
         covariates = "X", outcome_surrogate = "S_Y", mc_draws = 0)
}
