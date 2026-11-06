# MissingMediatorAnalysis

Mediation analysis with a partially missing mediator and outcome.

`MissingMediatorAnalysis` estimates the exposure → mediator → outcome indirect effect (product of
coefficients) when the mediator, and optionally the outcome, are missing at random given fully
observed synthetic surrogates (predictions of the missing variables from data available for
everyone). The proposed method, **SynMediation**, extends synthetic surrogate (SynSurr) analysis to
mediation and uses all subjects through three least-squares fits:

1. surrogates on the exposure and covariates (all *n* subjects);
2. mediator on surrogates, exposure and covariates (subjects with the mediator);
3. outcome on mediator, surrogates, exposure and covariates (subjects with both).

Inference uses the stacked influence functions of both paths, so the covariance between the two
paths is accounted for; the indirect effect has delta-method and Monte Carlo confidence intervals.
The outcome can be continuous (linear model, default) or binary (`outcome_family = "binomial"`,
logistic model; see below). The package depends only on base R.

## Installation

```r
# install.packages("remotes")
remotes::install_github("haoyu-yang001/MissingMediatorAnalysis")
```

or from a local copy: `R CMD INSTALL MissingMediatorAnalysis`.

## Functions

| Function | Purpose |
|---|---|
| `mma()` | Full mediation analysis from a data frame (SynMediation by default, `outcome_method = "SeqSynSurr"`). |
| `mma_scan()` | Many exposures × mediators (e.g. thousands of proteins), one row per pair; errors and warnings are recorded instead of stopping; optional parallel cores. |
| `mma_diagnostics()` | Checks of the working assumptions for one pathway (surrogate strength, linearity, constant variance, missingness model, comparison with complete case, nonlinear sensitivity fit). |
| `fit_synsurr()` | Exposure–mediator path alone (SynSurr or complete case). |
| `fit_outcome_model()` | Alternative mediator–outcome estimators (SurrAIPW, SynSurr, Combined, AIPWnew, AIPW, IPW, CC). |
| `sim_mma_data()` | Simulate data with a missing mediator (and outcome). |
| `print`, `summary`, `coef`, `vcov`, `as.data.frame` | Methods for `"mma"` objects. |

## Input

`mma()` takes one `data.frame` with one row per subject; columns are referred to by name.

| Role | Argument | Requirements |
|---|---|---|
| Exposure | `exposure` | numeric, complete; continuous or 0/1 |
| Mediator | `mediator` | numeric, **may contain `NA`** |
| Outcome | `outcome` | numeric, **may contain `NA`** if `outcome_surrogate` is given |
| Mediator surrogate | `surrogate` | numeric, complete; prediction of the mediator. May be a vector: the surrogate first, then fully observed auxiliary variables on which missingness may depend |
| Outcome surrogate | `outcome_surrogate` | numeric, complete; prediction of the outcome; required if the outcome has `NA` |
| Covariates | `covariates` | character vector of column names, complete |

Rows with missing exposure, surrogates or covariates are dropped with a message
(`na_action = "omit"`, default) or cause an error (`na_action = "fail"`). Impute small amounts of
missing surrogate inputs (e.g. by the median) before computing the surrogates.

## Output

`mma()` returns a list of class `"mma"`. Main components:

| Component | Content |
|---|---|
| `effects` | rows `alpha`, `beta`, `NIE`, `NDE`, `TE`, `PM`; columns `estimate`, `se`, `z`, `p`, `lower`, `upper` |
| `nie_ci` | NIE confidence intervals: delta method and Monte Carlo (recommended) |
| `tests` | tests of NIE = 0: Sobel, joint significance, Monte Carlo |
| `vcov`, `cor_ab` | joint covariance of the six effects; estimated corr(alpha, beta) |
| `mediator_model`, `outcome_model` | all coefficients of the two models with SE, z, p and CI |
| `n`, `n_obs`, `n_obs_Y`, `n_obs_both` | total; with mediator; with outcome; with both |
| `weight_diagnostics`, `missing_data` | weighted SynMediation (`weighting = "ipw"`): weight diagnostics and missingness models |

## Example

```r
library(MissingMediatorAnalysis)
set.seed(1)
dat <- sim_mma_data(n = 10000, miss_rate = 0.8, rho = 0.8, mechanism = "MAR",
                    rho_Y = 0.8, miss_rate_Y = 0.5, mechanism_Y = "MAR")
fit <- mma(dat, exposure = "A", mediator = "M", outcome = "Y",
           surrogate = "S", covariates = "X", outcome_surrogate = "S_Y")
summary(fit)
```

```
      estimate      se     z         p    lower   upper
alpha  0.12240 0.01603 7.633 2.297e-14 0.090950 0.15380
beta   0.13350 0.02130 6.267 3.686e-10 0.091730 0.17520
NIE    0.01633 0.00344 4.748 2.059e-06 0.009591 0.02308
...
```

Only 1,369 of the 10,000 simulated subjects have both the mediator and the outcome; all 10,000
contribute. Without `outcome_surrogate` (fully observed outcome) the same call runs
SynMediation with the mediator surrogate only.

### Weighted SynMediation (sensitivity analysis)

The closed form above needs, under missingness at random given the surrogates, that steps 2 and 3
are linear regressions. `weighting = "ipw"` removes this requirement: steps 2 and 3 are fitted by
weighted least squares with inverse probabilities of being observed (logistic missingness models
P(mediator observed | surrogates, exposure, covariates) in everyone and
P(outcome observed | mediator, surrogates, exposure, covariates) among subjects with the mediator),
and all formulas are otherwise unchanged. The estimator then targets the same associations as a
complete-case analysis under data missing completely at random, whatever the relation between the
surrogates and the measurements; its validity depends only on the missingness models, and with
equal weights it is the closed form. Expect larger standard errors.

```r
fit_w <- mma(dat, exposure = "A", mediator = "M", outcome = "Y", surrogate = "S",
             covariates = "X", outcome_surrogate = "S_Y", weighting = "ipw")
fit_w$weight_diagnostics   # largest weights, effective sample sizes, smallest probabilities,
                           # largest single-subject shares of Var(alpha) and Var(beta)
```

The missingness models use the `surrogate` columns (including auxiliary columns such as an
attendance score), `outcome_surrogate`, the exposure and the covariates (plus extra `ps_vars`); the
outcome is never used (ordering mediator -> outcome; subjects with the outcome but no mediator count as
missing the outcome). With a complete outcome only the mediator model is fitted. `ps_trim` bounds
each fitted probability; `basis = "poly2"`/`"ns3"` is supported. When one subject carries more than
20% of Var(alpha) or Var(beta) the fit warns; refit with `ps_trim` and without that subject.
`mma_scan(..., weighting = "ipw")` adds the weight diagnostics to each row.

### Scanning many exposures

```r
res <- mma_scan(data, exposures = G,          # n x (number of proteins) matrix
                mediators  = c("VAT", "ASAT"),
                surrogates = c("S_VAT", "S_ASAT"),
                outcome = "FVC", covariates = covs,
                outcome_surrogate = "S_FVC",
                mc_draws = 1e4, n_cores = 32, seed = 1)
mma_diagnostics(fit, data)
```

### Binary outcome

For a 0/1 outcome (e.g. COPD), set `outcome_family = "binomial"`. The outcome model is logistic,
so `beta`, `NIE`, `NDE` and `TE` are on the log odds ratio scale and `fit$effects_or` gives odds
ratios. Both paths target the full-cohort projections (the complete-case targets under MCAR) and
are estimated by a sequential augmented inverse probability weighted (AIPW) estimator
(`outcome_method = "SeqSynSurr"`, default). Its validity depends only on the logistic missingness
models; the surrogates and the working prediction models only affect precision.
`outcome_method = "IPW"` gives the weighted complete-case fit with the same weights and `"CC"` the
complete-case fit.

The missing-data ordering (`missing_order`) is chosen automatically and stated by `print()`:
with a complete outcome, the probability that the mediator is observed is modelled on (surrogates, exposure,
covariates, outcome); with a monotone pattern (outcome recorded only with the mediator) the order
is M -> Y; otherwise Y -> M, where whether the mediator is observed may depend on the recorded outcome.
`fit$binary_diagnostics` reports events, the largest weight, the effective sample size, the
smallest propensities and the largest single-subject share of Var(beta). With fewer than about 50
events among subjects with mediator and outcome, or when one subject carries more than 20% of
Var(beta), the fit warns; report sensitivity analyses with `ps_trim` and without that subject.

```r
fit <- mma(data, exposure = "protein", mediator = "MFI", outcome = "COPD",
           surrogate = c("S_MFI", "attendance_score"), covariates = covs,
           outcome_surrogate = "S_COPD", outcome_family = "binomial", mc_draws = 1e4)
fit$effects_or
fit$binary_diagnostics
```

`mma_scan()` passes `outcome_family` and the other arguments to `mma()`; its output includes
`outcome_family` and `events`. `mma_diagnostics()` is available for continuous outcomes only.

## Reference

Yang H, Wang X, Lin X. MRI fat compartments beyond BMI link the plasma proteome to lung function:
surrogate-assisted mediation analysis in UK Biobank. Manuscript in preparation.

## License

GPL (>= 3)
