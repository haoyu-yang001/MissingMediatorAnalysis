\name{mma}
\alias{mma}
\alias{print.mma}
\alias{summary.mma}
\alias{coef.mma}
\alias{vcov.mma}
\title{Mediation analysis with a partially missing mediator and outcome}
\description{
Product-of-coefficients mediation analysis of
\eqn{M = \alpha_0 + \alpha_A A + X^\top\alpha_X + \epsilon_M} and
\eqn{Y = \beta_0 + \beta_A A + \beta_M M + X^\top\beta_X + \epsilon_Y} when the
mediator, and optionally the outcome, are partially missing and synthetic
surrogates are observed for everyone. By default both paths are estimated by
SynMediation (\code{outcome_method = "SeqSynSurr"}), the maximum
likelihood estimator of the factorized Gaussian model
\eqn{f(S \mid Z) f(M \mid S, Z) f(Y \mid M, S, Z)}, with \eqn{S} the surrogate
of \eqn{M} (and of \eqn{Y}, if given) and \eqn{Z = (1, A, X)}. The function
reports \eqn{NIE = \alpha_A\beta_M}, \eqn{NDE = \beta_A}, \eqn{TE = NDE + NIE}
and the proportion mediated. Standard errors use the stacked influence
functions of both paths, so they include \eqn{Cov(\hat\alpha_A, \hat\beta_M)}.
Two tests of \eqn{H_0: NIE = 0} are given: the delta-method (Sobel-type) Wald
test and the joint-significance test.
}
\usage{
mma(data, exposure, mediator, outcome, surrogate, covariates = NULL,
    ps_vars = NULL, mediator_method = c("SynSurr", "CC"),
    outcome_method = c("SeqSynSurr", "Combined", "SurrAIPW", "SynSurr",
                       "AIPWnew", "AIPW", "IPW", "CC"),
    conf_level = 0.95, ps_trim = NULL, mc_draws = 1e5, outcome_surrogate = NULL,
    na_action = c("omit", "fail"), keep_IF = TRUE, verbose = TRUE,
    basis = c("linear", "poly2", "ns3"), basis_outcome = FALSE,
    outcome_family = c("gaussian", "binomial"),
    missing_order = c("auto", "YM", "MY"), firth = TRUE, gh_nodes = 15)
}
\arguments{
\item{data}{a data frame.}
\item{exposure, mediator, outcome, surrogate}{column names. The mediator may
contain \code{NA}, and so may the outcome if \code{outcome_surrogate} is given.
For \code{"SeqSynSurr"}, \code{surrogate} may name several columns: the
synthetic surrogate of the mediator first, then auxiliary variables (fully
observed, not in \code{covariates}) on which missingness may depend; they are
integrated out like the surrogates.}
\item{covariates}{character vector of covariate column names; factor and
character columns are expanded into indicators. Covariates that are constant
within the subset with an observed mediator (or with both observed) are not
identified there; their coefficients are reported as \code{NA}, which does not
affect \eqn{\alpha_A} or \eqn{\beta_M}.}
\item{ps_vars}{variables (besides an intercept) for the propensity model and the working model of \eqn{E[M \mid W]}; defaults to \code{c(outcome, exposure, covariates, surrogate)}, or \code{c(exposure, covariates, surrogate, outcome_surrogate)} when the outcome has missing values. Not used by \code{"SeqSynSurr"} or \code{"SynSurr"} for a continuous outcome. For a binary outcome the missingness models always contain the surrogates, \code{outcome_surrogate}, the exposure and the covariates (and the outcome as set by \code{missing_order}); further fully observed variables named in \code{ps_vars} are added to them like auxiliary surrogate columns.}
\item{mediator_method}{\code{"SynSurr"} or complete-case \code{"CC"}.}
\item{outcome_method}{\code{"SeqSynSurr"} (default; requires
\code{mediator_method = "SynSurr"}, whose fit it replaces by SynSurr with all
surrogates). Alternatives with a complete outcome: \code{"Combined"} (optimal
combination of AIPWnew and CC; requires that missingness does not depend on the
outcome given the regressors), \code{"AIPWnew"}, \code{"AIPW"}, \code{"IPW"} or
\code{"CC"}. Alternatives with a missing outcome: \code{"SurrAIPW"}
(surrogate-augmented AIPW for monotone missingness), \code{"SynSurr"} (SynSurr
for the outcome among subjects with an observed mediator; requires that the
mediator's missingness does not depend on the outcome given the regressors) or
\code{"CC"}. With \code{outcome_family = "binomial"}: \code{"SeqSynSurr"}
(sequential AIPW, default), \code{"IPW"} (weighted complete case) or
\code{"CC"}; the first two estimate both paths from the same inverse-probability
weights (see Details).}
\item{conf_level}{confidence level.}
\item{ps_trim}{optional lower bound for fitted propensities. For a binary
outcome it is applied to each missingness model separately (a sensitivity
analysis for extreme weights; the influence functions ignore it).}
\item{mc_draws}{number of Monte Carlo draws for the NIE confidence interval;
\code{0} skips it.}
\item{na_action}{\code{"omit"} (default) drops, with a message, rows with
missing values in the exposure, covariates, surrogates or \code{ps_vars};
\code{"fail"} stops instead.}
\item{keep_IF}{keep the influence-function matrices in \code{fits}; set to
\code{FALSE} to save memory when fitting many pathways.}
\item{verbose}{print messages about dropped rows.}
\item{basis_outcome}{with \code{basis = "poly2"} or \code{"ns3"}, also use the
basis of the surrogates in step 3 (sensitivity analysis; less efficient).}
\item{basis}{for \code{"SeqSynSurr"}: \code{"linear"} (default) or a
nonlinear model for \eqn{E(M \mid S, A, X)} in step 2, \code{"poly2"} (adds
the square of each surrogate) or \code{"ns3"} (natural spline with 3 df per
surrogate); step 3 stays linear. The paths are then obtained by averaging the
fitted moments over all subjects. Use it as a sensitivity analysis when
\code{\link{mma_diagnostics}} flags nonlinearity in step 2. Only the average
residual variance of step 2 enters the estimator, estimated among subjects with
\eqn{M}; if this variance depends on the surrogates, check the
\code{heteroskedasticity} part of \code{\link{mma_diagnostics}}. With a strongly
heteroskedastic step 2 the nonlinear bases are not necessarily more reliable
than the linear fit.}
\item{outcome_family}{\code{"gaussian"} (default, linear outcome model) or
\code{"binomial"} for a 0/1 outcome (logistic outcome model; see Details).
With \code{"binomial"}, \code{outcome_method} must be \code{"SeqSynSurr"},
\code{"IPW"} or \code{"CC"} and \code{basis} must be \code{"linear"}.}
\item{missing_order}{for a binary outcome with \code{"SeqSynSurr"} or
\code{"IPW"}: \code{"auto"} (default), \code{"YM"} or \code{"MY"}, the
ordering of the sequential missing-at-random assumption (see Details).
\code{"auto"} uses \code{"MY"} when every subject with the outcome also has the
mediator and the outcome is partially missing, and \code{"YM"} otherwise. The
ordering and the missingness models are stated by \code{print()}.}
\item{firth}{for a binary outcome, use Firth's penalized likelihood in the
logistic fits among subjects with mediator and outcome: the weighted fit of
\code{"IPW"} (which also starts the AIPW iterations), the working model of the
outcome for \code{missing_order = "MY"}, and \code{"CC"}. The final AIPW
equation is solved over all subjects and is not penalized.}
\item{gh_nodes}{number of Gauss-Hermite points used to integrate over the
working model of a missing mediator (binary outcome).}
\item{outcome_surrogate}{optional column name of a synthetic surrogate of the
outcome (e.g. an earlier measurement), observed for everyone. With a missing
outcome it is required by \code{"SurrAIPW"} and \code{"SynSurr"}; for
\code{"SeqSynSurr"} it is optional and weakens the MAR assumption for the
outcome (missingness may then depend on it).}
}
\details{
\code{"SeqSynSurr"} fits (1) each surrogate on \eqn{Z} in all subjects, (2)
\eqn{M} on \eqn{(S, Z)} among subjects with an observed mediator and (3)
\eqn{Y} on \eqn{(M, S, Z)} among subjects with both observed, and recovers
\eqn{\alpha} and the regression of \eqn{Y} on \eqn{(Z, M)} by integrating over
\eqn{S}. It is valid when the mediator is missing at random given
\eqn{(A, X, S)} and the outcome given \eqn{(A, X, S, M)}, so missingness may
depend on the surrogates, and when the conditional means in steps 2 and 3 are
linear with constant residual variance in step 2. Subjects with an outcome but
no mediator are treated as missing the outcome. With no missing data it
reproduces least squares exactly.

\strong{Binary outcome.} With \code{outcome_family = "binomial"} the outcome
model is \eqn{logit P(Y = 1 \mid A, X, M) = \beta_0 + \beta_A A + \beta_M M +
X^\top\beta_X}. Both paths are defined as full-cohort projections (least
squares of \eqn{M} on \eqn{Z = (1, A, X)}; logistic projection of \eqn{Y} on
\eqn{(Z, M)}), which are also the complete-case targets when data are missing
completely at random, so \eqn{\beta_M} is a log odds ratio per unit of the
mediator and \eqn{NIE = \alpha_A\beta_M} is on the log odds ratio scale. The
default (\code{outcome_method = "SeqSynSurr"}) is a sequential augmented
inverse probability weighted (AIPW) estimator. Its consistency rests only on
the missingness models: the surrogates and the working
prediction models for \eqn{M} (normal linear in the surrogates and \eqn{Z},
clipped to the observed range of \eqn{M}) and for \eqn{Y} (logistic) affect
precision, never the target. Both paths use the same missingness models:
\eqn{\alpha} is the least-squares fit of \eqn{\hat m + w(M - \hat m)} on
\eqn{Z}, and \eqn{(\beta_Z, \beta_M)} solves, over all subjects, the
logistic score equation in which predicted scores are corrected by
inverse-probability-weighted residuals of the subjects with the mediator and
outcome. \code{outcome_method = "IPW"} is the weighted complete-case fit with
the same weights (Firth-penalized weighted logistic regression) and
\code{"CC"} the unweighted complete-case fit.

With \eqn{S} the \code{surrogate} columns, \code{outcome_surrogate} and any
extra \code{ps_vars}, and \eqn{R_M}, \eqn{R_Y} the indicators that the mediator
and the outcome are observed, the logistic missingness models are, by
\code{missing_order}:
\describe{
  \item{\code{"YM"}}{(default when the outcome is complete, or when some
    subjects have the outcome but not the mediator)
    \eqn{P(R_Y = 1 \mid S, Z)} and \eqn{P(R_M = 1 \mid S, Z, R_Y, R_Y Y)}, so
    whether the mediator is observed may depend on the recorded outcome. With a
    complete outcome this is \eqn{P(R_M = 1 \mid S, Z, Y)}.}
  \item{\code{"MY"}}{(monotone; default when every subject with the outcome
    also has the mediator) \eqn{P(R_M = 1 \mid S, Z)} and
    \eqn{P(R_Y = 1 \mid R_M = 1, M, S, Z)}; subjects with the outcome but no
    mediator are treated as missing the outcome.}
}
The influence functions include the working-model fits and the estimation of
the missingness models, and are stacked with those of \eqn{\alpha} as for a
continuous outcome. Simulations show reliable inference from about 50 events
among subjects with mediator and outcome. The fit stores
\code{binary_diagnostics}: \code{events_MY} (events among subjects with
mediator and outcome), \code{events_Y} (among all with the outcome),
\code{max_weight}, \code{ess} (effective sample size of the weights),
\code{min_pi_M}, \code{min_pi_Y} (smallest fitted probabilities of being observed
among the observed), \code{converged}, \code{eq_norm}, and
\code{max_if_share_beta}, the largest single-subject share of the estimated
variance of \eqn{\hat\beta_M}, with the row name of that subject
(\code{most_influential}). It warns, also in \code{print()}, when there are
fewer than 50 events or one subject carries more than 20\% of that variance;
then report sensitivity analyses with \code{ps_trim} and without the most
influential subject. The function stops when there are no events and warns when
there are fewer than 5 per parameter. The output adds \code{effects_or}, the
odds ratios of \eqn{\beta_M}, NIE, NDE and TE. With no missing data the fit
reproduces \code{lm} and \code{glm(..., family = binomial)} exactly.

The tests of \eqn{H_0: NIE = 0} are the delta-method (Sobel-type) Wald
test, the joint-significance test \eqn{\max(p_\alpha, p_\beta)} and the
p-value matching the Monte Carlo interval,
\eqn{2\min\{P(\alpha^*\beta^* < 0), P(\alpha^*\beta^* > 0)\}} under the
bivariate normal approximation, computed by numerical integration (so it is
not limited by the number of draws). The function stops when the number of
subjects with the mediator, or with both mediator and outcome, does not exceed
the number of parameters, and warns when it is below 10 per parameter.

Besides the delta-method (Wald) interval, a Monte Carlo confidence interval for
the NIE is reported. It draws \eqn{(\alpha_A^*, \beta_M^*)} from the bivariate
normal approximation to \eqn{(\hat\alpha_A, \hat\beta_M)}, including their
estimated covariance, and takes quantiles of \eqn{\alpha_A^*\beta_M^*}, so the
interval reflects the skewness of a product of two estimates. It uses the random
number generator; call \code{set.seed()} beforehand for reproducible limits.
}
\value{An object of class \code{"mma"} with components \code{effects},
\code{tests}, \code{nie_ci} (delta-method and Monte Carlo confidence intervals
for the NIE), \code{conf_level}, \code{vcov}, \code{cor_ab},
\code{mediator_model}, \code{outcome_model}, \code{fits}, \code{n},
\code{n_obs} (observed mediator), \code{n_obs_Y} (observed outcome),
\code{n_obs_both}, \code{outcome_family} and, for a binary outcome,
\code{events} (events among subjects with mediator and outcome),
\code{effects_or}, \code{binary_diagnostics} and \code{missing_data} (the
ordering used and a description of the missingness models).}
\examples{
dat <- sim_mma_data(n = 2000, miss_rate = 0.5, rho = 0.8)
fit <- mma(dat, exposure = "A", mediator = "M", outcome = "Y",
           surrogate = "S", covariates = "X")
summary(fit)
# Mediator and outcome both partially missing
dat2 <- sim_mma_data(n = 2000, miss_rate = 0.75, rho = 0.8,
                     rho_Y = 0.8, miss_rate_Y = 0.5)
mma(dat2, exposure = "A", mediator = "M", outcome = "Y", surrogate = "S",
    covariates = "X", outcome_surrogate = "S_Y")
# Binary outcome (sequential AIPW)
dat2$Yb <- ifelse(is.na(dat2$Y), NA, as.numeric(dat2$Y > quantile(dat2$Y, 0.9, na.rm = TRUE)))
mma(dat2, exposure = "A", mediator = "M", outcome = "Yb", surrogate = "S",
    covariates = "X", outcome_surrogate = "S_Y", outcome_family = "binomial",
    mc_draws = 1e4)
}
