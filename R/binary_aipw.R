# Binary outcome: target-preserving sequential AIPW (and IPW) estimators.
#
# Targets (full-cohort projections; the complete-case targets under MCAR):
#   alpha: least squares of M on Z = (1, A, X);
#   theta = (beta_Z, beta_M): logistic projection of Y on x = (Z, M), the solution of
#          E[x {Y - expit(x'theta)}] = 0;  NIE = alpha_A beta_M (log odds ratio scale).
# Validity rests only on the missingness models; the surrogates and the
# working prediction models (for M and for Y) affect precision, never the target.
#
# Two orderings of sequential missingness at random (S = surrogates and auxiliaries):
#   "YM" (.fit_binary_ym): R_Y depends on (S, Z); R_M on (S, Z, R_Y, R_Y * Y). With a
#        complete outcome this is the "P2" pattern, P(R_M = 1 | S, Z, Y).
#        theta solves sum_i [h0 + wY (h1 - h0) + wYM (u - h1)] = 0 with
#        wY = R_Y / pi_Y, wYM = R_Y R_M / (pi_Y pi_M), u = x {Y - expit(x'theta)},
#        h1 = E~[u | S, Z, Y] (Gauss-Hermite over M ~ N(m(S, Z, 1, Y), s2)) and
#        h0 = p_Y h1(Y = 1) + (1 - p_Y) h1(Y = 0) (p_Y: weighted ML logistic of Y on
#        (S, Z) among R_Y = 1). alpha = OLS of m + wM (M - m) on Z, wM = R_M / pi_M.
#   "MY" (.fit_binary_my): monotone, Y used only with M ("P1"). R_M depends on (S, Z);
#        R_Y, among R_M = 1, on (M, S, Z). theta solves
#        sum_i [h0 + w1 (h1 - h0) + w2 (u - h1)] = 0, w1 = R_M / pi_M,
#        w2 = R_M R_Y / (pi_M lambda_Y), h1 = x {mu3(M, S, Z) - expit(x'theta)} (mu3:
#        Firth logistic of Y on (M, S, Z) among R_M R_Y = 1), h0 = Gauss-Hermite mean of h1
#        over M ~ N(m(S, Z), s2). alpha = OLS of m + w1 (M - m) on Z.
# "IPW" is the weighted complete-case fit with the same weights (WLS of M on Z; Firth
# weighted logistic of Y on x).
# Influence functions (n scale, theta_hat - theta = n^-1 sum IF_i) include the
# working-model fits (finite-sample correction; their expected contribution is zero
# under correct missingness models) and are residualised on the scores of the missingness
# models. The working prediction of M is clipped to the observed range of M.
# Ported from binary_outcome/robust/R/synmed_robust.R (method "AIPW"/"IPW") and
# binary_outcome/robust/ymiss/R/synmed_ymiss.R.

# (Weighted) logistic regression on the rows with w > 0, optionally with Firth's
# penalty (weighted hat matrix), returning n-scale influence functions of the
# unpenalized weighted score (the penalty is O(1)). Columns aliased within the
# fitted rows get coefficient and influence function 0.
.logit_if <- function(X, y, w = rep(1, nrow(X)), firth = TRUE, maxit = 1000, tol = 1e-8) {
  n <- nrow(X)
  if (any(w < 0)) stop("negative weights")
  use <- w > 0
  Xu <- X[use, , drop = FALSE]; yu <- y[use]; wu <- w[use]
  if (sum(yu) == 0 || sum(yu) == length(yu))
    stop(sprintf("the outcome has no %s among the %d subjects with mediator and outcome observed",
                 if (sum(yu) == 0) "events" else "non-events", length(yu)))
  qx <- qr(Xu, tol = 1e-7)
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  Xk <- Xu[, keep, drop = FALSE]
  b <- numeric(ncol(Xk)); conv <- FALSE
  for (it in seq_len(maxit)) {
    p <- plogis(drop(Xk %*% b)); W <- wu * p * (1 - p)
    Iinv <- tryCatch(solve(crossprod(Xk * W, Xk)), error = function(e)
      stop("logistic fit failed (separation?)", if (!firth) "; use firth = TRUE", call. = FALSE))
    U <- wu * (yu - p)
    if (firth) U <- U + W * rowSums((Xk %*% Iinv) * Xk) * (0.5 - p)
    step <- drop(Iinv %*% crossprod(Xk, U))
    mx <- max(abs(step)); if (mx > 2) step <- step * 2 / mx
    b <- b + step
    if (mx < tol) { conv <- TRUE; break }
  }
  if (!conv) warning(sprintf("logistic fit did not converge in %d iterations", maxit), call. = FALSE)
  p <- plogis(drop(Xk %*% b))
  J <- crossprod(Xk * (wu * p * (1 - p)), Xk) / n
  coef <- numeric(ncol(X)); coef[keep] <- b; names(coef) <- colnames(X)
  IF <- matrix(0, n, ncol(X)); colnames(IF) <- colnames(X)
  IF[use, keep] <- (Xk * (wu * (yu - p))) %*% solve(J)
  list(coef = coef, IF = IF, converged = conv,
       aliased = colnames(X)[setdiff(seq_len(ncol(X)), keep)],
       events = sum(yu), n = sum(use))
}

# Probabilists' Gauss-Hermite rule: sum_q w_q f(z_q) approximates E f(N(0, 1)).
.gh_rule <- function(Q) {
  i <- seq_len(Q - 1)
  J <- matrix(0, Q, Q); J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1, ]^2)
}

# Central-difference Jacobian of fun at par (columns = parameters).
.numjac <- function(fun, par, eps = 1e-6) {
  h <- eps * (1 + abs(par))
  sapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, h[j])
    (fun(par + e) - fun(par - e)) / (2 * h[j])
  })
}

# Logistic missingness model P(R = 1 | V) fitted on the rows `use`: fitted values
# (1 elsewhere) and the score matrix (0 outside `use`).
.bin_propensity <- function(V, R, use = rep(TRUE, length(R))) {
  n <- length(R); pr <- rep(1, n); sc <- matrix(0, n, ncol(V))
  f <- suppressWarnings(glm.fit(V[use, , drop = FALSE], as.numeric(R[use]), family = binomial()))
  if (!f$converged) warning("missingness (propensity) model did not converge", call. = FALSE)
  pr[use] <- f$fitted.values
  sc[use, ] <- (R[use] - pr[use]) * V[use, , drop = FALSE]
  list(p = pr, score = sc)
}

# Residualise influence functions on the missingness-model scores (estimation of the
# propensities); no-op without a fitted missingness model.
.bin_res <- function(IF, score) if (is.null(score)) IF else .residualize(IF, score)

# Keep the intercept and the non-constant columns.
.bin_dropconst <- function(V)
  V[, c(TRUE, apply(V[, -1, drop = FALSE], 2, function(v) diff(range(v)) > 0)), drop = FALSE]

# Ordering M -> Y (monotone, "P1"). Returns alpha, theta with influence functions
# and weight diagnostics.
.fit_binary_my <- function(Y, M, S, Z, method = c("AIPW", "IPW"), firth = TRUE, ps_trim = NULL,
                           Q = 15, clip_m = TRUE, nuisance_if = TRUE) {
  method <- match.arg(method)
  n <- nrow(Z); p <- ncol(Z)
  R1 <- !is.na(M); R2 <- R1 & !is.na(Y)
  Mz <- ifelse(R1, M, 0); Yz <- ifelse(R2, Y, 0)

  ## missingness models and weights
  score <- NULL; pi1 <- rep(1, n); lam <- rep(1, n)
  if (!all(R1)) {
    ps1 <- .bin_propensity(cbind(S, Z), R1)
    pi1 <- ps1$p; score <- ps1$score
  }
  lam_fit <- any(R1 & !R2)
  if (lam_fit) {
    ps2 <- .bin_propensity(cbind(M = Mz, S, Z), R2, use = R1)
    lam <- ps2$p; score <- cbind(score, ps2$score)
  }
  # Trimming is applied to each model's probability before combining; the influence
  # functions ignore it.
  if (!is.null(ps_trim)) { pi1 <- pmax(pi1, ps_trim); lam <- pmax(lam, ps_trim) }
  w1 <- R1 / pi1; w2 <- R2 / (pi1 * lam)
  x <- cbind(Z, M = Mz)

  ## weighted complete-case fits (IPW; also the Newton start)
  fa <- .wls_if(Z, Mz, w = w1)
  fb <- .logit_if(x, Yz, w2, firth = firth, maxit = 500, tol = 1e-9)
  diag <- list(n_M = sum(R1), n_Y = sum(!is.na(Y)), n_MY = sum(R2), events_MY = sum(Yz[R2]),
               max_weight = max(w2), ess = sum(w2)^2 / sum(w2^2),
               min_pi_M = if (!all(R1)) min(pi1[R1]) else NA_real_,
               min_pi_Y = if (lam_fit) min(lam[R2]) else NA_real_,
               converged = NA, eq_norm = NA_real_)
  if (method == "IPW")
    return(list(alpha = fa$coef, IF_a = .bin_res(fa$IF, score), theta = fb$coef,
                IF_b = .bin_res(fb$IF, score), diag = diag, aliased = fb$aliased))

  ## working models
  gh <- .gh_rule(Q)
  W0m <- cbind(S, Z)
  fm <- .wls_if(W0m, Mz, w = as.numeric(R1))
  Mrange <- if (clip_m) range(Mz[R1]) else c(-Inf, Inf)
  clipm <- function(v) pmin(pmax(v, Mrange[1]), Mrange[2])
  m <- clipm(drop(W0m %*% fm$coef)); s2 <- sum(fm$resid[R1]^2) / sum(R1)
  Mq <- m + sqrt(s2) * outer(rep(1, n), gh$z)
  X3 <- cbind(M = Mz, S, Z)
  f3 <- .logit_if(X3, Yz, w = as.numeric(R2), firth = firth)
  cSZ <- drop(cbind(S, Z) %*% f3$coef[-1])
  Yq <- plogis(f3$coef[1] * Mq + cSZ); Yhat1 <- plogis(f3$coef[1] * Mz + cSZ)
  dmu <- function(e) { q <- plogis(e); q * (1 - q) }
  h0_fun <- function(th) {
    Rq <- Yq - plogis(drop(Z %*% th[1:p]) + th[p + 1] * Mq)
    cbind(Z * drop(Rq %*% gh$w), M = drop((Mq * Rq) %*% gh$w))
  }

  ## AIPW: sum_i [(1 - w1) h0 + w1 x (Yhat1 - mu) + w2 x (Y - Yhat1)] = 0, Newton.
  th <- fb$coef; conv <- FALSE
  for (it in 1:100) {
    eq <- drop(Z %*% th[1:p]) + th[p + 1] * Mq; V <- dmu(eq)
    cvec <- 1 - w1
    a <- drop(V %*% gh$w); b <- drop((Mq * V) %*% gh$w); c2 <- drop((Mq^2 * V) %*% gh$w)
    zb <- crossprod(Z, cvec * b)
    H <- rbind(cbind(crossprod(Z * (cvec * a), Z), zb), c(zb, sum(cvec * c2)))
    H <- H + crossprod(x * (w1 * dmu(drop(x %*% th))), x)
    r1 <- w1 * (Yhat1 - plogis(drop(x %*% th))) + w2 * (Yz - Yhat1)
    psi <- cvec * h0_fun(th) + x * r1
    step <- drop(solve(H, colSums(psi)))
    mx <- max(abs(step)); if (mx > 2) step <- step * 2 / mx
    th <- th + step
    if (mx < 1e-10) { conv <- TRUE; break }
  }
  if (!conv) warning("AIPW estimating equation did not converge (max step ", signif(mx, 3), ")",
                     call. = FALSE)
  r1 <- w1 * (Yhat1 - plogis(drop(x %*% th))) + w2 * (Yz - Yhat1)
  psi_b <- (1 - w1) * h0_fun(th) + x * r1
  diag$converged <- conv
  diag$eq_norm <- sqrt(sum(colSums(psi_b)^2)) / n
  fa_aipw <- .wls_if(Z, m + w1 * (Mz - m))
  if (nuisance_if) {
    Ja <- crossprod(Z) / n
    xi <- c(fm$coef, s2, f3$coef)
    IF_xi <- cbind(fm$IF, R1 * (fm$resid^2 - s2) * n / sum(R1), f3$IF)
    q0 <- length(fm$coef)
    psi_xi <- function(xi) {
      mm <- clipm(drop(W0m %*% xi[1:q0])); M_q <- mm + sqrt(xi[q0 + 1]) * outer(rep(1, n), gh$z)
      om <- xi[q0 + 1 + seq_along(f3$coef)]; cc <- drop(cbind(S, Z) %*% om[-1])
      Y_q <- plogis(om[1] * M_q + cc); Yh <- plogis(om[1] * Mz + cc)
      Rq <- Y_q - plogis(drop(Z %*% th[1:p]) + th[p + 1] * M_q)
      hh <- cbind(Z * drop(Rq %*% gh$w), drop((M_q * Rq) %*% gh$w))
      rr <- w1 * (Yh - plogis(drop(x %*% th))) + w2 * (Yz - Yh)
      colMeans((1 - w1) * hh + x * rr)
    }
    psi_b <- psi_b + IF_xi %*% t(.numjac(psi_xi, xi, eps = 1e-5))
    Da_xi <- crossprod(Z * (1 - w1), W0m) / n
    fa_aipw$IF <- fa_aipw$IF + (fm$IF %*% t(Da_xi)) %*% t(solve(Ja))
  }
  list(alpha = fa_aipw$coef, IF_a = .bin_res(fa_aipw$IF, score), theta = th,
       IF_b = .bin_res(psi_b %*% t(solve(H / n)), score), diag = diag, aliased = character())
}

# Ordering Y -> M (general; "P2" when the outcome is complete).
.fit_binary_ym <- function(Y, M, S, Z, method = c("AIPW", "IPW"), firth = TRUE, ps_trim = NULL,
                           Q = 15, clip_m = TRUE) {
  method <- match.arg(method)
  n <- nrow(Z); p <- ncol(Z)
  RY <- !is.na(Y); RM <- !is.na(M)
  Ycomp <- all(RY); Mcomp <- all(RM)
  Yz <- ifelse(RY, Y, 0); Mz <- ifelse(RM, M, 0)
  x <- cbind(Z, M = Mz)

  ## missingness models and weights
  score <- NULL; piY <- rep(1, n); piM <- rep(1, n)
  if (!Ycomp) { ps <- .bin_propensity(cbind(S, Z), RY); piY <- ps$p; score <- ps$score }
  VM <- .bin_dropconst(cbind(Z, S, R_Y = as.numeric(RY), RY_Y = Yz))   # (S, Z, R_Y, R_Y Y)
  if (!Mcomp) { ps <- .bin_propensity(VM, RM); piM <- ps$p; score <- cbind(score, ps$score) }
  if (!is.null(ps_trim)) { piY <- pmax(piY, ps_trim); piM <- pmax(piM, ps_trim) }
  wY <- RY / piY; wM <- RM / piM; wYM <- RY * RM / (piY * piM)

  fipw <- .logit_if(x, Yz, wYM, firth = firth, maxit = 500, tol = 1e-9)
  diag <- list(n_M = sum(RM), n_Y = sum(RY), n_MY = sum(RY & RM), events_MY = sum(Yz[RY & RM]),
               max_weight = max(wYM), ess = sum(wYM)^2 / sum(wYM^2),
               min_pi_M = if (!Mcomp) min(piM[RM]) else NA_real_,
               min_pi_Y = if (!Ycomp) min(piY[RY]) else NA_real_,
               converged = NA, eq_norm = NA_real_)
  if (method == "IPW") {
    fwa <- .wls_if(Z, Mz, w = wM)
    return(list(alpha = fwa$coef, IF_a = .bin_res(fwa$IF, score), theta = fipw$coef,
                IF_b = .bin_res(fipw$IF, score), diag = diag, aliased = fipw$aliased))
  }

  ## working models
  gh <- .gh_rule(Q); one <- rep(1, n)
  if (!Ycomp) {
    VY <- cbind(Z, S, if (Mcomp) cbind(M = Mz))
    fp <- .logit_if(VY, Yz, wY, firth = FALSE, maxit = 500, tol = 1e-9)
    pY_fun <- function(b) plogis(drop(VY %*% b))
  }
  if (!Mcomp) {
    fm <- .wls_if(VM, Mz, w = as.numeric(RM))
    Mrange <- if (clip_m) range(Mz[RM]) else c(-Inf, Inf)
    clipm <- function(v) pmin(pmax(v, Mrange[1]), Mrange[2])
    s2 <- sum(fm$resid[RM]^2) / sum(RM)
    VMat <- function(ry, yy) cbind(Z, S, R_Y = ry, RY_Y = yy)[, colnames(VM), drop = FALSE]
    V1y <- VMat(1, Yz); V11 <- VMat(1, 1); V10 <- VMat(1, 0)
  }
  # Node sets (M values, Y values, node weights) for h1 (level 1) and h0 (level 0).
  nodes <- function(xi) {
    if (!Mcomp) {
      sd <- sqrt(xi$s2); zq <- outer(one, gh$z); wq <- outer(one, gh$w)
      L1 <- list(M = clipm(drop(V1y %*% xi$m)) + sd * zq, Y = matrix(Yz, n, Q), W = wq)
      if (Ycomp) return(list(L1 = L1, L0 = L1))
      pr <- pY_fun(xi$p)
      L0 <- list(M = cbind(clipm(drop(V11 %*% xi$m)) + sd * zq, clipm(drop(V10 %*% xi$m)) + sd * zq),
                 Y = cbind(matrix(1, n, Q), matrix(0, n, Q)), W = cbind(pr * wq, (1 - pr) * wq))
    } else {
      pr <- pY_fun(xi$p)
      L1 <- list(M = cbind(Mz), Y = cbind(Yz), W = cbind(one))
      L0 <- list(M = cbind(Mz, Mz), Y = cbind(one, 0 * one), W = cbind(pr, 1 - pr))
    }
    list(L1 = L1, L0 = L0)
  }
  hfun <- function(L, th) {
    Rq <- L$Y - plogis(drop(Z %*% th[1:p]) + th[p + 1] * L$M)
    cbind(Z * rowSums(L$W * Rq), M = rowSums(L$W * L$M * Rq))
  }
  hjac <- function(L, th, cc) {            # sum_i cc_i * (- d h_i / d theta)
    mu <- plogis(drop(Z %*% th[1:p]) + th[p + 1] * L$M); V <- L$W * mu * (1 - mu)
    a <- rowSums(V); b <- rowSums(V * L$M); c2 <- rowSums(V * L$M^2)
    zb <- crossprod(Z, cc * b)
    rbind(cbind(crossprod(Z * (cc * a), Z), zb), c(zb, sum(cc * c2)))
  }
  xi0 <- list(m = if (!Mcomp) fm$coef, s2 = if (!Mcomp) s2, p = if (!Ycomp) fp$coef)
  c0 <- 1 - wY; c1 <- wY - wYM; c2 <- wYM
  psi_fun <- function(th, nd) {
    mu <- plogis(drop(x %*% th))
    c0 * hfun(nd$L0, th) + c1 * hfun(nd$L1, th) + c2 * x * (Yz - mu)
  }

  ## beta: Newton from the weighted complete-case (Firth) fit
  nd <- nodes(xi0); th <- fipw$coef; conv <- FALSE
  for (it in 1:100) {
    mu <- plogis(drop(x %*% th))
    H <- hjac(nd$L0, th, c0) + hjac(nd$L1, th, c1) + crossprod(x * (c2 * mu * (1 - mu)), x)
    step <- drop(solve(H, colSums(psi_fun(th, nd))))
    mx <- max(abs(step)); if (mx > 2) step <- step * 2 / mx
    th <- th + step
    if (mx < 1e-10) { conv <- TRUE; break }
  }
  if (!conv) warning("AIPW estimating equation did not converge (max step ", signif(mx, 3), ")",
                     call. = FALSE)
  mu <- plogis(drop(x %*% th))
  H <- hjac(nd$L0, th, c0) + hjac(nd$L1, th, c1) + crossprod(x * (c2 * mu * (1 - mu)), x)
  psi_b <- psi_fun(th, nd)
  diag$converged <- conv
  diag$eq_norm <- sqrt(sum(colSums(psi_b)^2)) / n

  ## alpha: OLS of m + wM (M - m) on Z
  if (!Mcomp) {
    m_own <- clipm(drop(VM %*% fm$coef))
    fa <- .wls_if(Z, m_own + wM * (Mz - m_own))
  } else fa <- .wls_if(Z, Mz)

  ## working-model estimation in the influence functions
  pack <- function(v) {
    i <- 0; take <- function(k) { r <- v[i + seq_len(k)]; i <<- i + k; r }
    list(m = if (!Mcomp) take(length(fm$coef)), s2 = if (!Mcomp) take(1),
         p = if (!Ycomp) take(length(fp$coef)))
  }
  xv <- c(xi0$m, xi0$s2, xi0$p)
  IF_xi <- cbind(if (!Mcomp) cbind(fm$IF, RM * (fm$resid^2 - s2) * n / sum(RM)), if (!Ycomp) fp$IF)
  Dxi <- .numjac(function(v) colMeans(psi_fun(th, nodes(pack(v)))), xv, eps = 1e-5)
  IF_b <- (psi_b + IF_xi %*% t(Dxi)) %*% t(solve(H / n))
  if (!Mcomp) {
    Da <- crossprod(Z * (1 - wM), VM) / n
    fa$IF <- fa$IF + (fm$IF %*% t(Da)) %*% t(solve(crossprod(Z) / n))
  }
  list(alpha = fa$coef, IF_a = .bin_res(fa$IF, score), theta = th,
       IF_b = .bin_res(IF_b, score), diag = diag, aliased = character())
}

# Binary-outcome fit used by mma(): resolves the missing-data ordering, calls the
# engine and returns both paths in the format of the other outcome fits.
.fit_binary <- function(Y, M, S, Z, mediator_name, method = c("AIPW", "IPW"),
                        missing_order = c("auto", "YM", "MY"), firth = TRUE, ps_trim = NULL,
                        Q = 15, row_ids = NULL) {
  method <- match.arg(method); missing_order <- match.arg(missing_order)
  S <- as.matrix(S)
  n <- nrow(Z); p <- ncol(Z)
  RM <- !is.na(M); RY <- !is.na(Y)
  sv <- paste(colnames(S), collapse = ", ")
  order <- if (all(RM) && all(RY)) "none" else if (missing_order != "auto") missing_order else
    if (all(RY) || any(RY & !RM)) "YM" else "MY"
  if (order == "none") {
    fa <- .wls_if(Z, M)
    fb <- .logit_if(cbind(Z, M = M), Y, firth = method == "IPW" && firth)
    res <- list(alpha = fa$coef, IF_a = fa$IF, theta = fb$coef, IF_b = fb$IF,
                diag = list(n_M = n, n_Y = n, n_MY = n, events_MY = sum(Y), max_weight = 1,
                            ess = n, min_pi_M = NA_real_, min_pi_Y = NA_real_,
                            converged = fb$converged, eq_norm = NA_real_),
                aliased = fb$aliased)
    assumption <- "no missing mediator or outcome: full-data least squares and logistic regression"
  } else if (order == "YM") {
    res <- .fit_binary_ym(Y, M, S, Z, method, firth, ps_trim, Q)
    assumption <- if (all(RY))
      sprintf("outcome complete; P(mediator observed) logistic in (%s, exposure, covariates, outcome)", sv) else
      if (all(RM)) sprintf("mediator complete; P(outcome observed) logistic in (%s, exposure, covariates)", sv) else
      sprintf("sequential MAR, order Y -> M: P(outcome observed) logistic in (%s, exposure, covariates); P(mediator observed) logistic in (%s, exposure, covariates, R_Y, R_Y * outcome)", sv, sv)
  } else {
    res <- .fit_binary_my(Y, M, S, Z, method, firth, ps_trim, Q)
    assumption <- sprintf("sequential MAR, order M -> Y: P(mediator observed) logistic in (%s, exposure, covariates)%s",
                          sv, if (any(RM & !RY))
                            sprintf("; P(outcome observed | mediator observed) logistic in (mediator, %s, exposure, covariates)", sv) else "")
    if (any(RY & !RM)) assumption <- paste0(assumption, "; subjects with the outcome but no mediator are treated as missing the outcome")
  }
  xn <- c(colnames(Z), mediator_name)
  bad <- intersect(res$aliased, xn[c(2, p + 1)])
  if (length(bad)) stop("not identified among subjects with mediator and outcome: ", paste(bad, collapse = ", "))
  theta <- setNames(unname(res$theta), xn); IF_b <- res$IF_b; colnames(IF_b) <- xn
  alpha <- setNames(unname(res$alpha), colnames(Z)); IF_a <- res$IF_a; colnames(IF_a) <- colnames(Z)
  ib <- IF_b[, p + 1]
  d <- res$diag
  d$events_Y <- sum(Y, na.rm = TRUE)
  d$max_if_share_beta <- max(ib^2) / sum(ib^2)
  d$most_influential <- if (is.null(row_ids)) which.max(ib^2) else row_ids[which.max(ib^2)]
  label <- if (order == "none") "full data" else
    sprintf("%s (%s)", if (method == "AIPW") "sequential AIPW" else "IPW", gsub("(.)(.)", "\\1 -> \\2", order))
  if (method == "IPW" && firth && order != "none") label <- sub(")", ", Firth)", label, fixed = TRUE)
  list(coef = theta, vcov = .vcov_if(IF_b), IF = IF_b, method = label,
       mediator = list(coef = alpha, IF = IF_a, vcov = .vcov_if(IF_a), method = label),
       diagnostics = d, missing_order = order, assumption = assumption,
       aliased = intersect(res$aliased, xn[-c(1, 2, p + 1)]), n_obs_both = sum(RM & RY))
}

# Complete-case logistic regression of Y on (Z, M) among subjects with both.
.fit_logit_cc <- function(Y, M, Zt, firth = TRUE) {
  R2 <- as.numeric(!is.na(M) & !is.na(Y))
  X <- Zt; X[R2 == 0, ncol(X)] <- 0
  f <- .logit_if(X, ifelse(R2 == 1, Y, 0), w = R2, firth = firth)
  list(coef = f$coef, vcov = .vcov_if(f$IF), IF = f$IF,
       method = if (firth) "CC logistic (Firth)" else "CC logistic",
       steps = list(events = f$events), aliased = f$aliased, n_obs_both = sum(R2))
}
