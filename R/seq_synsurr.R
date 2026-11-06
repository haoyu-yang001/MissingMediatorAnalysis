# Sequential SynSurr for the mediator-outcome path when both the mediator and
# the outcome are partially missing.
#
# With surrogates S (n x k, e.g. of the mediator and of the outcome) observed
# for everyone and a monotone pattern (M observed: R1; M and Y observed: R2),
# the Gaussian likelihood factorizes as f(S | Z) f(M | S, Z) f(Y | M, S, Z).
# Under MAR (R1 depends on (Z, S), R2 given R1 = 1 on (Z, S, M)) each factor is
# fitted by least squares on the subjects whose response is observed:
#   1. S on Z, all subjects:            kappa (p x k), Sigma_S
#   2. M on (S, Z) among R1 = 1:        delta, gamma, sigma2_M
#   3. Y on (M, S, Z) among R2 = 1:     omega_M, omega_S, omega_Z
# and the outcome regression Y ~ Z + M is recovered by integrating out S given
# (Z, M): E(S | M, Z) = kappa' Z + nu (M - Z' alpha), with
# nu = Sigma_S delta / (delta' Sigma_S delta + sigma2_M), so
#   alpha = gamma + kappa delta,  beta_M = omega_M + nu' omega_S,
#   beta_Z = omega_Z + (kappa - alpha nu') omega_S.
# Subjects with Y but no M are treated as missing Y (monotone pattern).
# Influence functions of the least-squares fits and the residual (co)variances
# are stacked and mapped by the delta method (numerical Jacobian).
#
# Weighted SynMediation: w1 (step 2, zero without M) and w2 (step 3, zero
# without M and Y) replace the indicators R1 and R2; steps 2 and 3 are then
# weighted least squares and sigma2_M is the weighted mean squared step-2
# residual. With w1 = R1 / pi_M and w2 = R1 R2 / (pi_M lambda_Y) the fits
# estimate the full-cohort projections and the estimator is target-preserving
# (see .fit_seq_ipw, which also accounts for the estimated weights). Weights
# proportional to R1 and R2 give exactly the unweighted estimator.
.fit_seq_synsurr <- function(Y, M, S, Z, Zt, basis = "linear", basis_outcome = FALSE,
                             w1 = NULL, w2 = NULL) {
  S <- as.matrix(S)
  n <- nrow(Z); p <- ncol(Z); k <- ncol(S)
  if (is.null(colnames(S)) || any(colnames(S) == "")) colnames(S) <- paste0("S", seq_len(k))
  if (basis != "linear") return(.fit_seq_basis(Y, M, S, Z, Zt, basis, basis_outcome, w1, w2))
  R1 <- as.numeric(!is.na(M))
  R2 <- R1 * !is.na(Y)
  if (is.null(w1)) w1 <- R1
  if (is.null(w2)) w2 <- R2
  Mz <- ifelse(R1 == 1, M, 0)   # zero-filled M only enters rows with weight 0

  f1 <- lapply(seq_len(k), function(j) .wls_if(Z, S[, j]))
  kappa <- sapply(f1, `[[`, "coef")
  E_S <- sapply(f1, `[[`, "resid")
  Sigma_S <- crossprod(E_S) / n
  f2 <- .wls_if(cbind(S, Z), M, w = w1)
  sigma2_M <- sum(w1 * f2$resid^2) / sum(w1)
  f3 <- .wls_if(cbind(M = Mz, S, Z), Y, w = w2)
  # The exposure, the surrogates and the mediator must be identified in every
  # fit; covariates aliased within a subset (e.g. a centre with no subject
  # observed) get coefficient 0, which leaves alpha_A and beta unaffected.
  key <- c(colnames(S), colnames(Z)[2], "M")
  bad <- intersect(c(unlist(lapply(f1, `[[`, "aliased")), f2$aliased, f3$aliased), key)
  if (length(bad))
    stop("not identified (collinear within the observed subset): ", paste(bad, collapse = ", "))
  aliased <- unique(c(f2$aliased, f3$aliased))

  # Parameter vector: kappa, lower triangle of Sigma_S, (delta, gamma), sigma2_M,
  # (omega_M, omega_S, omega_Z).
  lt <- lower.tri(Sigma_S, diag = TRUE)
  par <- c(kappa, Sigma_S[lt], f2$coef, sigma2_M, f3$coef)
  map <- function(par) {
    i <- 0
    take <- function(m) { v <- par[i + seq_len(m)]; i <<- i + m; v }
    kap <- matrix(take(p * k), p, k)
    Sig <- matrix(0, k, k); Sig[lt] <- take(sum(lt)); Sig <- Sig + t(Sig) - diag(diag(Sig), k)
    dg <- take(k + p); delta <- dg[1:k]; gamma <- dg[-(1:k)]
    s2 <- take(1)
    om <- take(1 + k + p); omM <- om[1]; omS <- om[1 + 1:k]; omZ <- om[-(1:(1 + k))]
    alpha <- gamma + drop(kap %*% delta)
    nu <- drop(Sig %*% delta) / drop(crossprod(delta, Sig %*% delta) + s2)
    bS <- sum(omS * nu)
    c(alpha, omZ + drop(kap %*% omS) - alpha * bS, omM + bS)
  }

  # Influence functions of the parameters, on the scale of all n subjects.
  IF_Sig <- do.call(cbind, lapply(which(lt), function(ix) {
    j <- row(Sigma_S)[ix]; l <- col(Sigma_S)[ix]
    E_S[, j] * E_S[, l] - Sigma_S[j, l]
  }))
  IF_s2 <- w1 * (f2$resid^2 - sigma2_M) * n / sum(w1)
  IF_par <- cbind(do.call(cbind, lapply(f1, `[[`, "IF")), IF_Sig, f2$IF, IF_s2, f3$IF)

  est <- map(par)
  delta <- f2$coef[seq_len(k)]
  nu <- drop(Sigma_S %*% delta) / drop(crossprod(delta, Sigma_S %*% delta) + sigma2_M)
  names(nu) <- colnames(S)
  # Partial R^2 of M on the surrogates given Z, among subjects with M
  # (weighted by w1 in weighted SynMediation).
  rz <- .wls_if(Z, M, w = w1)$resid
  partial_r2 <- 1 - sigma2_M / .wmean(rz^2, w1)
  h <- 1e-5 * ifelse(par == 0, 1, abs(par))   # relative step
  Jac <- sapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, h[j])
    (map(par + e) - map(par - e)) / (2 * h[j])
  })
  IF <- IF_par %*% t(Jac)
  ia <- seq_len(p); ib <- p + seq_len(p + 1)
  coef_b <- setNames(est[ib], colnames(Zt))
  IF_b <- IF[, ib, drop = FALSE]; colnames(IF_b) <- colnames(Zt)
  IF_a <- IF[, ia, drop = FALSE]; colnames(IF_a) <- colnames(Z)
  list(coef = coef_b, vcov = .vcov_if(IF_b), IF = IF_b, method = "SeqSynSurr",
       mediator = list(coef = setNames(est[ia], colnames(Z)), IF = IF_a,
                       vcov = .vcov_if(IF_a), method = "SynSurr"),
       steps = list(kappa = kappa, Sigma_S = Sigma_S, delta_gamma = f2$coef,
                    sigma2_M = sigma2_M, omega = f3$coef, nu = nu,
                    beta_M_direct = f3$coef[["M"]],
                    beta_M_surrogate = sum(nu * f3$coef[1 + seq_len(k)]),
                    partial_r2 = partial_r2),
       aliased = aliased, n_obs_both = sum(R2))
}

# Weighted mean (the plain mean over the rows with w = 1 for 0/1 weights).
.wmean <- function(v, w) if (all(w == 0 | w == 1)) mean(v[w == 1]) else sum(w * v) / sum(w)

# Sequential SynSurr with a nonlinear mediator model in step 2,
# E(M | S, Z) = B(S)'b + Z'g with B a quadratic ("poly2") or natural spline
# ("ns3", 3 df) basis of each surrogate; step 3 is linear in S, or uses the
# same basis with basis_outcome = TRUE (h = B(S)'w_B + Z'w_Z). The two
# paths are then the moments
#   alpha = E[ZZ']^{-1} E[Z m],  (beta_Z, beta_M) = A^{-1} c,
#   A = [E[ZZ'], E[Zm]; E[mZ'], E[m^2] + sigma2_M],
#   c = (E[Z(w_M m + h)], E[w_M (m^2 + sigma2_M) + m h]),
# with m the fitted step-2 mean, h = S'w_S + Z'w_Z, and expectations averaged
# over all n subjects. With B(S) = S this equals the linear version. Only the
# average residual variance enters, so a heteroskedastic step 2 biases beta
# only through (mean over all - mean over observed) * (w_M - beta_M).
# Influence function: averaging with the working models fixed, plus the
# estimation of the step-2 and step-3 fits and of sigma2_M (numerical Jacobian).
.fit_seq_basis <- function(Y, M, S, Z, Zt, basis = c("poly2", "ns3"), basis_outcome = FALSE,
                           w1 = NULL, w2 = NULL) {
  basis <- match.arg(basis)
  n <- nrow(Z); p <- ncol(Z); k <- ncol(S)
  B <- switch(basis,
    poly2 = { S2 <- S^2; colnames(S2) <- paste0(colnames(S), "^2"); cbind(S, S2) },
    ns3 = do.call(cbind, lapply(seq_len(k), function(j) {
      b <- unclass(splines::ns(S[, j], df = 3))[, 1:3, drop = FALSE]
      colnames(b) <- paste0(colnames(S)[j], "_ns", 1:3); b })))
  R1 <- as.numeric(!is.na(M)); R2 <- R1 * !is.na(Y)
  if (is.null(w1)) w1 <- R1
  if (is.null(w2)) w2 <- R2
  Mz <- ifelse(R1 == 1, M, 0)

  X2 <- cbind(B, Z); f2 <- .wls_if(X2, M, w = w1)
  sigma2_M <- sum(w1 * f2$resid^2) / sum(w1)
  IF_s2 <- w1 * (f2$resid^2 - sigma2_M) * n / sum(w1)
  X3 <- cbind(M = Mz, if (basis_outcome) B else S, Z); f3 <- .wls_if(X3, Y, w = w2)
  key <- c(colnames(Z)[2], "M", if (basis == "poly2") colnames(S))
  bad <- intersect(c(f2$aliased, f3$aliased), key)
  if (length(bad))
    stop("not identified (collinear within the observed subset): ", paste(bad, collapse = ", "))
  q2 <- ncol(X2); q3 <- ncol(X3)
  par <- c(f2$coef, sigma2_M, f3$coef)
  Ezz <- crossprod(Z) / n

  parts <- function(par) {
    m <- drop(X2 %*% par[1:q2]); s2 <- par[q2 + 1]; om <- par[q2 + 1 + 1:q3]
    list(m = m, s2 = s2, wM = om[1], h = drop(X3[, -1, drop = FALSE] %*% om[-1]))
  }
  moments <- function(q) {
    Ezm <- drop(crossprod(Z, q$m)) / n
    list(Ezm = Ezm, A = rbind(cbind(Ezz, Ezm), c(Ezm, mean(q$m^2) + q$s2)),
         c = c(drop(crossprod(Z, q$wM * q$m + q$h)) / n,
               mean(q$wM * (q$m^2 + q$s2) + q$m * q$h)))
  }
  est_fun <- function(par) {
    mo <- moments(parts(par))
    c(solve(Ezz, mo$Ezm), solve(mo$A, mo$c))
  }
  est <- est_fun(par)
  ia <- seq_len(p); ib <- p + seq_len(p + 1)
  alpha <- est[ia]; theta <- est[ib]

  q <- parts(par); mo <- moments(q)
  lz <- drop(Z %*% theta[1:p]); tM <- theta[p + 1]
  g <- cbind(Z * (q$wM * q$m + q$h - lz - q$m * tM),
             q$wM * (q$m^2 + q$s2) + q$m * q$h - q$m * lz - (q$m^2 + q$s2) * tM)
  IF_fixed <- cbind((Z * (q$m - drop(Z %*% alpha))) %*% t(solve(Ezz)), g %*% t(solve(mo$A)))
  IF_par <- cbind(f2$IF, IF_s2, f3$IF)
  h <- 1e-5 * ifelse(par == 0, 1, abs(par))
  Jac <- sapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, h[j])
    (est_fun(par + e) - est_fun(par - e)) / (2 * h[j])
  })
  IF <- IF_fixed + IF_par %*% t(Jac)
  IF_a <- IF[, ia, drop = FALSE]; colnames(IF_a) <- colnames(Z)
  IF_b <- IF[, ib, drop = FALSE]; colnames(IF_b) <- colnames(Zt)
  rz <- .wls_if(Z, M, w = w1)$resid
  list(coef = setNames(theta, colnames(Zt)), vcov = .vcov_if(IF_b), IF = IF_b,
       method = paste0("SeqSynSurr(", basis, if (basis_outcome) ", both steps", ")"),
       mediator = list(coef = setNames(alpha, colnames(Z)), IF = IF_a,
                       vcov = .vcov_if(IF_a), method = paste0("SynSurr(", basis, ")")),
       steps = list(mean_model = f2$coef, sigma2_M = sigma2_M, omega = f3$coef,
                    beta_M_direct = f3$coef[["M"]],
                    beta_M_surrogate = theta[p + 1] - f3$coef[["M"]],
                    partial_r2 = 1 - sigma2_M / .wmean(rz^2, w1)),
       aliased = unique(c(f2$aliased, f3$aliased)), n_obs_both = sum(R2))
}

# Weighted SynMediation (continuous outcome, ordering M -> Y). Missingness
# models (logistic, maximum likelihood):
#   pi_M = P(R_M = 1 | S, Z) in all n subjects (if M is partially missing);
#   lambda_Y = P(R_Y = 1 | R_M = 1, M, S, Z) among subjects with M (if some of
#   them lack Y). Subjects with Y but no M are treated as missing Y.
# Steps 2 and 3 of .fit_seq_synsurr are weighted by w1 = R_M / pi_M and
# w2 = R_M R_Y / (pi_M lambda_Y) (w2 = w1 with a complete outcome); step 1 and
# the closed-form formulas are unchanged. Every moment is then a full-cohort
# moment estimated by inverse probability weighting, so the estimator targets
# the full-cohort projections (the complete-case targets under MCAR) whatever
# the relation between S and (M, Y); validity depends only on the missingness
# models. It equals the sequential AIPW estimator with these weighted working
# fits (all augmentation terms vanish). The influence functions are those of
# the weighted closed form, residualised on the scores of the missingness
# models (estimated weights); ps_trim bounds each fitted probability from
# below (the influence functions ignore the trimming).
.fit_seq_ipw <- function(Y, M, S, Z, Zt, basis = "linear", basis_outcome = FALSE,
                         ps_trim = NULL, row_ids = NULL) {
  S <- as.matrix(S)
  n <- nrow(Z); p <- ncol(Z)
  R1 <- !is.na(M); R2 <- R1 & !is.na(Y)
  Mz <- ifelse(R1, M, 0)
  score <- NULL; pi1 <- rep(1, n); lam <- rep(1, n); conv <- TRUE
  pi_fit <- !all(R1); lam_fit <- any(R1 & !R2)
  if (pi_fit) {
    ps1 <- .bin_propensity(cbind(S, Z), R1)
    pi1 <- ps1$p; score <- ps1$score; conv <- ps1$converged
  }
  if (lam_fit) {
    ps2 <- .bin_propensity(cbind(M = Mz, S, Z), R2, use = R1)
    lam <- ps2$p; score <- cbind(score, ps2$score); conv <- conv && ps2$converged
  }
  if (!is.null(ps_trim)) { pi1 <- pmax(pi1, ps_trim); lam <- pmax(lam, ps_trim) }
  w1 <- R1 / pi1; w2 <- R2 / (pi1 * lam)
  fit <- .fit_seq_synsurr(Y, M, S, Z, Zt, basis, basis_outcome, w1 = w1, w2 = w2)
  if (!is.null(score)) {
    fit$IF[] <- .residualize(fit$IF, score)
    fit$mediator$IF[] <- .residualize(fit$mediator$IF, score)
    fit$vcov <- .vcov_if(fit$IF); fit$mediator$vcov <- .vcov_if(fit$mediator$IF)
  }
  order <- if (!pi_fit && !lam_fit) "none" else "MY"
  sv <- paste(colnames(S), collapse = ", ")
  assumption <- if (order == "none")
    "no missing mediator or outcome: weights are 1 (full-data least squares)" else
    paste0(if (pi_fit) sprintf("P(mediator observed) logistic in (%s, exposure, covariates)", sv) else
             "mediator complete",
           if (lam_fit) sprintf("; P(outcome observed | mediator observed) logistic in (mediator, %s, exposure, covariates)", sv) else
             if (all(!is.na(Y))) "; outcome complete" else
               "; every subject with the mediator has the outcome",
           if (any(!is.na(Y) & !R1) && !all(!is.na(Y)))
             "; subjects with the outcome but no mediator are treated as missing the outcome" else "")
  share <- function(v) max(v^2) / sum(v^2)
  ia <- fit$mediator$IF[, 2]; ib <- fit$IF[, p + 1]
  rid <- function(i) if (is.null(row_ids)) i else row_ids[i]
  fit$diagnostics <- list(
    n_M = sum(R1), n_Y = sum(!is.na(Y)), n_MY = sum(R2),
    max_weight_M = max(w1), ess_M = sum(w1)^2 / sum(w1^2),
    max_weight = max(w2), ess = sum(w2)^2 / sum(w2^2),
    min_pi_M = if (pi_fit) min(pi1[R1]) else NA_real_,
    min_pi_Y = if (lam_fit) min(lam[R2]) else NA_real_,
    max_if_share_alpha = share(ia), most_influential_alpha = rid(which.max(ia^2)),
    max_if_share_beta = share(ib), most_influential = rid(which.max(ib^2)),
    converged = conv)
  fit$method <- paste0(sub("^SeqSynSurr", "SeqSynSurr-IPW", fit$method),
                       if (order == "MY") " (M -> Y)" else "")
  fit$mediator$method <- paste0(sub("^SynSurr", "SynSurr-IPW", fit$mediator$method))
  fit$missing_order <- order
  fit$assumption <- assumption
  fit$weights <- list(w1 = w1, w2 = w2)
  fit
}
