# Internal helpers shared by the mediator and outcome fits.
#
# Every estimator in the package is an M-estimator, and each fit returns an
# n x p matrix of influence functions IF such that
#   theta_hat - theta = n^{-1} sum_i IF_i + o_p(n^{-1/2}).
# Variances, and covariances across the two mediation paths, are then
# crossprod(IF) / n^2.

# Weighted least squares solving sum_i w_i x_i (y_i - x_i' b) = 0.
# Rows with w_i = 0 do not contribute; their y (and x) may be NA.
# Columns that are aliased among the contributing rows (e.g. an indicator that
# is zero for all of them) are detected by a pivoted QR decomposition; their
# coefficients and influence functions are set to 0 and their names returned
# in `aliased`. Columns are pivoted out from the right, so earlier columns are
# kept in preference to later ones.
.wls_if <- function(X, y, w = rep(1, nrow(X)), tol = 1e-7) {
  n <- nrow(X)
  use <- w != 0
  Xu <- X[use, , drop = FALSE]
  qx <- qr(Xu * sqrt(abs(w[use])), tol = tol)
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  Xk <- Xu[, keep, drop = FALSE]
  J <- crossprod(Xk * w[use], Xk) / n
  Jinv <- solve(J)
  b <- drop(Jinv %*% crossprod(Xk * w[use], y[use]) / n)
  coef <- numeric(ncol(X))
  coef[keep] <- b
  names(coef) <- colnames(X)
  res <- numeric(n)
  res[use] <- y[use] - drop(Xk %*% b)
  IF <- matrix(0, n, ncol(X))
  IF[use, keep] <- (Xk * (w[use] * res[use])) %*% Jinv
  colnames(IF) <- colnames(X)
  Jfull <- matrix(0, ncol(X), ncol(X))
  Jfull[keep, keep] <- J
  drop_idx <- setdiff(seq_len(ncol(X)), keep)
  aliased <- if (is.null(colnames(X))) drop_idx else colnames(X)[drop_idx]
  list(coef = coef, IF = IF, J = Jfull, resid = res, aliased = aliased)
}

.vcov_if <- function(IF) {
  V <- crossprod(IF) / nrow(IF)^2
  dimnames(V) <- list(colnames(IF), colnames(IF))
  V
}

# Moore-Penrose inverse; used for the (K + p) x (K + p) matrix B2 of the
# AIPW_new projection, which can be close to singular.
.pinv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  keep <- s$d > tol * s$d[1]
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# Logistic propensity model P(R = 1 | W), fitted by maximum likelihood.
.fit_propensity <- function(W, R, trim = NULL) {
  fit <- glm.fit(W, R, family = binomial())
  if (!fit$converged) warning("propensity model did not converge")
  omega <- fit$fitted.values
  if (!is.null(trim)) omega <- pmin(pmax(omega, trim), 1 - trim)
  list(omega = omega, eta = fit$coefficients,
       score = (R - omega) * W)
}

# Remove from psi its least-squares projection on the propensity score.
# This is the adjustment for estimating eta by maximum likelihood.
.residualize <- function(psi, score) {
  psi - score %*% .pinv(crossprod(score)) %*% crossprod(score, psi)
}

# Minimum-variance combination of two consistent estimators of the same
# parameter: theta_1 - Gamma (theta_1 - theta_2), with
# Gamma = Cov(theta_1, D) Var(D)^- and D = theta_1 - theta_2, all estimated from
# the influence functions. Gamma is treated as fixed in the returned IF.
.combine_if <- function(coef1, IF1, coef2, IF2) {
  D <- IF1 - IF2
  Gamma <- crossprod(IF1, D) %*% .pinv(crossprod(D))
  dimnames(Gamma) <- list(names(coef1), names(coef1))
  coef <- drop(coef1 - Gamma %*% (coef1 - coef2))
  names(coef) <- names(coef1)
  IF <- IF1 - D %*% t(Gamma)
  colnames(IF) <- colnames(IF1)
  list(coef = coef, IF = IF, Gamma = Gamma)
}

# Monte Carlo confidence interval for a * b: quantiles of the product of draws
# from the bivariate normal approximation N(est, V) to the two estimates.
.mc_product_ci <- function(est, V, draws, conf_level = 0.95) {
  e <- eigen(V, symmetric = TRUE)
  root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  z <- matrix(rnorm(2 * draws), draws, 2) %*% root
  ab <- (est[1] + z[, 1]) * (est[2] + z[, 2])
  q <- (1 - conf_level) / 2
  quantile(ab, c(q, 1 - q), names = FALSE)
}

.coef_table <- function(est, se, conf_level = 0.95) {
  est[!is.na(se) & se == 0] <- NA        # aliased (unidentified) coefficients
  se[is.na(est)] <- NA
  z <- est / se
  q <- qnorm(1 - (1 - conf_level) / 2)
  data.frame(estimate = est, se = se, z = z,
             p = 2 * pnorm(-abs(z)),
             lower = est - q * se, upper = est + q * se,
             row.names = names(est))
}

# P-value matching the Monte Carlo interval for a * b: 2 min{P(a* b* < 0),
# P(a* b* > 0)} under (a*, b*) ~ N(est, V). Writing a* = a + s_a z,
# b* | z ~ N(b + r s_b z, s_b^2 (1 - r^2)), the probability is a one-dimensional
# integral over z, split at the sign change of a*.
.product_p <- function(est, V) {
  sa <- sqrt(V[1, 1]); sb <- sqrt(V[2, 2])
  r <- max(min(V[1, 2] / (sa * sb), 1 - 1e-12), -1 + 1e-12)
  s <- sb * sqrt(1 - r^2)
  f <- function(z) {
    a <- est[1] + sa * z
    mu <- est[2] + r * sb * z
    dnorm(z) * ifelse(a > 0, pnorm(-mu / s), pnorm(mu / s))
  }
  z0 <- -est[1] / sa
  P <- stats::integrate(f, -Inf, z0, rel.tol = 1e-10, abs.tol = 0)$value +
    stats::integrate(f, z0, Inf, rel.tol = 1e-10, abs.tol = 0)$value
  P <- min(max(P, 0), 1)
  2 * min(P, 1 - P)
}

# Drop rows with missing values in variables that must be complete.
.drop_incomplete <- function(data, vars, na_action, verbose = TRUE) {
  bad <- !stats::complete.cases(data[, vars, drop = FALSE])
  if (!any(bad)) return(data)
  if (na_action == "fail")
    stop(sprintf("%d rows have missing values in %s", sum(bad),
                 paste(vars[colSums(is.na(data[, vars, drop = FALSE])) > 0], collapse = ", ")))
  if (verbose) message(sprintf("dropping %d of %d rows with missing values in %s", sum(bad),
                               nrow(data), paste(vars[colSums(is.na(data[, vars, drop = FALSE])) > 0],
                                                 collapse = ", ")))
  data[!bad, , drop = FALSE]
}

# Design matrix (intercept first) of the named variables; factor and character
# columns are expanded into indicators (treatment contrasts). The first
# variable keeps its name, so the exposure column is always named as given.
.design <- function(data, vars) {
  if (!length(vars)) return(cbind("(Intercept)" = rep(1, nrow(data))))
  f <- stats::reformulate(sprintf("`%s`", vars))
  X <- stats::model.matrix(f, data)
  colnames(X) <- gsub("`", "", colnames(X))
  if (is.numeric(data[[vars[1]]])) colnames(X)[2] <- vars[1]
  X
}

# Sample-size checks for the subsets used by each fit.
.check_sizes <- function(M, Y, p, k, outcome_method) {
  n_M <- sum(!is.na(M))
  n_both <- sum(!is.na(M) & !is.na(Y))
  if (n_M <= p + k)
    stop(sprintf("only %d subjects with an observed mediator for %d parameters", n_M, p + k))
  if (n_both <= p + k + 1)
    stop(sprintf("only %d subjects with mediator and outcome observed for %d parameters",
                 n_both, p + k + 1))
  if (n_both < 10 * (p + k + 1))
    warning(sprintf("%d subjects with mediator and outcome observed for %d parameters (< 10 per parameter); standard errors may be too small",
                    n_both, p + k + 1), call. = FALSE)
  invisible(NULL)
}
