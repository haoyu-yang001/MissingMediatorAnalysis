#' Diagnostics for SynMediation
#'
#' Checks the working assumptions of the default estimator for one pathway:
#' \describe{
#'   \item{surrogate}{partial \eqn{R^2} of the mediator on the surrogates given
#'     \eqn{Z}, and \eqn{R^2} of the mediator (outcome) on its own surrogate,
#'     among subjects with it observed.}
#'   \item{linearity}{robust (HC0) Wald tests of squared terms added to step 2
#'     (squares of the surrogates) and to step 3 (squares of the mediator and
#'     the surrogates).}
#'   \item{variance}{Breusch-Pagan test of constant residual variance in step 2
#'     (\eqn{nR^2} of the squared residuals on the surrogates and \eqn{Z}).}
#'   \item{missingness}{logistic models for observing the mediator given
#'     \eqn{(Z, S)} and for observing the outcome given \eqn{(Z, S, M)} among
#'     subjects with the mediator: Wald z of the surrogates (and mediator),
#'     area under the ROC curve and range of fitted probabilities.}
#'   \item{hausman}{test of SynMediation against complete-case analysis
#'     for \eqn{\alpha_A} and \eqn{\beta_M}, using the difference of their
#'     influence functions. Under MAR given the surrogates both are consistent
#'     only if missingness does not depend on the surrogates, so a significant
#'     difference is expected when it does; a large difference without
#'     surrogate-dependent missingness points to misspecification.}
#'   \item{heteroskedasticity}{the estimator uses the average step-2 residual
#'     variance among subjects with \eqn{M}; if the residual variance depends
#'     on \eqn{(S, Z)}, its average over all subjects (estimated by inverse
#'     probability weighting with the logistic model for observing \eqn{M})
#'     differs and \eqn{\hat\beta_M} changes by about
#'     \eqn{(\bar v_{all} - \bar v_{obs})(\omega_M - \beta_M) /
#'     Var(M \mid Z)}. Reported with its ratio to the standard error. The
#'     approximation assumes a correctly specified step-2 mean.}
#'   \item{nonlinearity}{\eqn{\alpha_A}, \eqn{\beta_M} and the NIE from
#'     the linear fit and with a quadratic step 2 (\code{basis = "poly2"}).}
#'   \item{surrogate_columns}{for each surrogate column (mediator surrogate,
#'     auxiliary columns, outcome surrogate): partial \eqn{R^2} on \eqn{Z},
#'     kurtosis of its step-1 residual and the share of the residual sum of
#'     squares from the 10 largest residuals, and \eqn{\hat\beta_M} and the ratio
#'     of standard errors when the column is left out. A large share (gross
#'     outliers in a covariate or surrogate) with a standard-error ratio above 1
#'     means the column adds noise rather than information; winsorize or drop it.}
#' }
#' All checks refer to the linear fit, also when \code{fit} used a nonlinear
#' \code{basis}.
#' @param fit an object of class \code{"mma"} fitted with
#'   \code{outcome_method = "SeqSynSurr"}.
#' @param data the data frame used for \code{fit}.
#' @return A list of data frames (class \code{"mma_diagnostics"}).
#' @export
mma_diagnostics <- function(fit, data) {
  v <- fit$vars
  if (identical(fit$outcome_family, "binomial")) stop("mma_diagnostics is implemented for a continuous outcome only")
  if (!startsWith(fit$methods[["outcome"]], "SeqSynSurr")) stop("diagnostics are for outcome_method = \"SeqSynSurr\"")
  data <- .drop_incomplete(as.data.frame(data),
                           unique(c(v$exposure, v$surrogate, v$covariates, v$outcome_surrogate)),
                           "omit", verbose = FALSE)
  M <- data[[v$mediator]]; Y <- data[[v$outcome]]
  Z <- .design(data, c(v$exposure, v$covariates))
  S <- cbind(as.matrix(data[, v$surrogate, drop = FALSE]),
             if (!is.null(v$outcome_surrogate)) data[[v$outcome_surrogate]])
  colnames(S) <- c(v$surrogate, v$outcome_surrogate)
  R1 <- !is.na(M); R2 <- R1 & !is.na(Y)
  Mz <- ifelse(R1, M, 0)
  r2 <- function(y, x, use) stats::cor(y[use], x[use])^2

  surrogate <- data.frame(
    quantity = c("partial R2 of M on surrogates | Z", paste("R2 of M on", v$surrogate[1]),
                 if (!is.null(v$outcome_surrogate)) paste("R2 of Y on", v$outcome_surrogate)),
    value = c(fit$fits$outcome$steps$partial_r2, r2(M, S[, 1], R1),
              if (!is.null(v$outcome_surrogate)) r2(Y, S[, ncol(S)], !is.na(Y))))

  # Robust Wald test of the last q columns of X in the fit of y on X (weights w).
  wald <- function(X, y, w, q) {
    f <- .wls_if(X, y, w)
    idx <- setdiff(ncol(X) - q + seq_len(q), match(f$aliased, colnames(X)))
    b <- f$coef[idx]
    V <- .vcov_if(f$IF)[idx, idx, drop = FALSE]
    stat <- drop(crossprod(b, solve(V, b)))
    c(df = length(idx), chisq = stat, p = stats::pchisq(stat, length(idx), lower.tail = FALSE))
  }
  S2 <- S^2; colnames(S2) <- paste0(colnames(S), "^2")
  lin2 <- wald(cbind(S, Z, S2), M, R1, ncol(S2))
  lin3 <- wald(cbind(M = Mz, S, Z, "M^2" = Mz^2, S2), Y, R2, 1 + ncol(S2))
  linearity <- data.frame(model = c("step 2: M ~ S + Z (+ S^2)", "step 3: Y ~ M + S + Z (+ M^2, S^2)"),
                          rbind(lin2, lin3), row.names = NULL)

  f2 <- .wls_if(cbind(S, Z), M, R1)
  e2 <- f2$resid[R1]^2
  X2 <- cbind(S, Z)[R1, , drop = FALSE]
  bp_r2 <- 1 - sum(stats::lm.fit(X2, e2)$residuals^2) / sum((e2 - mean(e2))^2)
  df_bp <- qr(X2)$rank - 1
  variance <- data.frame(test = "Breusch-Pagan, step 2", df = df_bp, chisq = sum(R1) * bp_r2,
                         p = stats::pchisq(sum(R1) * bp_r2, df_bp, lower.tail = FALSE))

  auc <- function(p, r) {
    rk <- rank(p); n1 <- as.numeric(sum(r)); n0 <- as.numeric(sum(!r))
    (sum(rk[r]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  logit <- function(X, r, terms) {
    g <- stats::glm.fit(X, as.numeric(r), family = stats::binomial())
    keep <- !is.na(g$coefficients)
    Xk <- X[, keep, drop = FALSE]
    V <- solve(crossprod(Xk * (g$weights), Xk))
    z <- g$coefficients[keep] / sqrt(diag(V))
    names(z) <- colnames(Xk)
    data.frame(term = terms, z = unname(z[terms]),
               auc = auc(g$fitted.values, r), min_prob = min(g$fitted.values),
               max_prob = max(g$fitted.values))
  }
  miss_M <- cbind(model = "observe M | Z, S", logit(cbind(Z, S), R1, colnames(S)))
  miss_Y <- if (any(R1 & !R2))
    cbind(model = "observe Y | Z, S, M (among M observed)",
          logit(cbind(Z, S, M = Mz)[R1, , drop = FALSE], R2[R1], c(colnames(S), "M"))) else NULL
  missingness <- rbind(miss_M, miss_Y)

  args <- list(data = data, exposure = v$exposure, mediator = v$mediator, outcome = v$outcome,
               surrogate = v$surrogate, covariates = v$covariates,
               outcome_surrogate = v$outcome_surrogate, mc_draws = 0, verbose = FALSE)
  ss <- suppressWarnings(do.call(mma, args))
  # complete case ignores the surrogates; pass only the mediator's own one
  cc <- suppressWarnings(do.call(mma, c(utils::modifyList(args, list(surrogate = v$surrogate[1])),
                                        mediator_method = "CC", outcome_method = "CC")))
  ifd <- function(f) cbind(alpha = f$fits$mediator$IF[, v$exposure], beta = f$fits$outcome$IF[, v$mediator])
  D <- ifd(ss) - ifd(cc)
  d <- c(alpha = ss$effects["alpha", "estimate"] - cc$effects["alpha", "estimate"],
         beta = ss$effects["beta", "estimate"] - cc$effects["beta", "estimate"])
  se_d <- sqrt(colSums(D^2)) / nrow(D)
  hausman <- data.frame(parameter = names(d), SeqSynSurr = c(ss$effects["alpha", "estimate"], ss$effects["beta", "estimate"]),
                        CC = c(cc$effects["alpha", "estimate"], cc$effects["beta", "estimate"]),
                        difference = d, se = se_d, z = d / se_d,
                        p = 2 * stats::pnorm(-abs(d / se_d)), row.names = NULL)

  X2a <- cbind(S, Z)
  pi1 <- stats::glm.fit(X2a, as.numeric(R1), family = stats::binomial())$fitted.values
  v_all <- sum(R1 * f2$resid^2 / pi1) / sum(R1 / pi1); v_obs <- mean(e2)
  m_all <- drop(X2a %*% f2$coef)
  var_mz <- mean(.wls_if(Z, m_all)$resid^2) + v_obs
  st <- ss$fits$outcome$steps
  b_lin <- ss$effects["beta", "estimate"]
  d_beta <- (v_all - v_obs) * (st$beta_M_direct - b_lin) / var_mz
  heteroskedasticity <- data.frame(mean_var_observed = v_obs, mean_var_all = v_all,
                                   implied_change_beta = d_beta,
                                   ratio_to_se = d_beta / ss$effects["beta", "se"])

  pq <- suppressWarnings(do.call(mma, c(args, basis = "poly2")))
  pick <- function(f) f$effects[c("alpha", "beta", "NIE"), "estimate"]
  nonlinearity <- data.frame(parameter = c("alpha", "beta", "NIE"), linear = pick(ss),
                             poly2 = pick(pq),
                             difference_in_se = (pick(pq) - pick(ss)) /
                               ss$effects[c("alpha", "beta", "NIE"), "se"],
                             row.names = NULL)

  # Per surrogate column: step-1 residual (on Z) diagnostics and the fit without
  # the column (the mediator's own surrogate is always kept).
  b0 <- ss$effects["beta", "se"]
  columns <- do.call(rbind, lapply(seq_len(ncol(S)), function(j) {
    e <- .wls_if(Z, S[, j])$resid
    e2s <- sort(e^2, decreasing = TRUE)
    f <- if (j == 1) NULL else {
      a <- args
      if (identical(colnames(S)[j], v$outcome_surrogate)) a["outcome_surrogate"] <- list(NULL)
      else a$surrogate <- setdiff(v$surrogate, colnames(S)[j])
      suppressWarnings(do.call(mma, a))
    }
    data.frame(column = colnames(S)[j], partial_R2_on_Z = 1 - mean(e^2) / stats::var(S[, j]),
               kurtosis_step1 = mean(e^4) / mean(e^2)^2,
               top10_share_step1 = sum(e2s[seq_len(min(10, length(e)))]) / sum(e2s),
               beta_without = if (is.null(f)) NA else f$effects["beta", "estimate"],
               se_ratio_with_vs_without = if (is.null(f)) NA else b0 / f$effects["beta", "se"])
  }))

  structure(list(surrogate = surrogate, linearity = linearity, variance = variance,
                 missingness = missingness, hausman = hausman,
                 heteroskedasticity = heteroskedasticity, nonlinearity = nonlinearity,
                 surrogate_columns = columns),
            class = "mma_diagnostics")
}

#' @export
print.mma_diagnostics <- function(x, digits = 3, ...) {
  for (nm in names(x)) {
    cat("\n", toupper(nm), "\n", sep = "")
    print(x[[nm]], digits = digits, row.names = FALSE)
  }
  invisible(x)
}
