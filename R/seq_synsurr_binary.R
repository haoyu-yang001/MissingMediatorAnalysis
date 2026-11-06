# SynMediation for a binary outcome.
#
# The mediator path is unchanged (steps 1-2 of .fit_seq_synsurr). The outcome
# model is logistic, logit P(Y = 1 | Z, M) = Z'beta + beta_M M, with
# (beta, beta_M) the full-cohort logistic projection, i.e. the solution of
# E[(Z, M)' {Y - expit(Z'beta + beta_M M)}] = 0.
#   3. logistic regression of Y on (M, S, Z) among R2 = 1 (Firth-penalized by
#      default: few events are a problem only in this fit).
#   4. "proj" (default): the full-cohort logistic regression of Y on (Z, M) is
#      solved after filling in what is missing. Subjects with M and Y contribute
#      their values; subjects with M but not Y contribute M and the step-3
#      probability; subjects without M contribute Q Gauss-Hermite points of the
#      step-2 normal model N(S'delta + Z'gamma, sigma2_M), each with the step-3
#      probability as a fractional outcome. Under MAR this estimating equation
#      has the same expectation as the full-data one.
#      "closed": the SynMediation closed form with a logistic step 3,
#      beta_M = omega_M + nu' omega_S, beta = omega_Z + (kappa - alpha nu') omega_S,
#      intercept + 0.5 omega_S' V omega_S, V = Var(S | Z, M). It is exact for a
#      log link and approximate for the logit; it overestimates |beta_M| as the
#      prevalence and the predictive strength of the surrogates increase.
# Influence functions: least-squares and logistic (maximum-likelihood form,
# evaluated at the Firth estimate, whose penalty is O(1)) influence functions of
# the steps, mapped by the delta method ("closed") or by the implicit-function
# theorem ("proj"), with numerical Jacobians.

# Logistic regression on the rows with w = 1 (w is a 0/1 indicator), optionally
# with Firth's penalty, returning n-scale influence functions. Columns aliased
# within the fitted rows get coefficient and influence function 0.
.logit_if <- function(X, y, w = rep(1, nrow(X)), firth = TRUE, maxit = 1000, tol = 1e-8) {
  n <- nrow(X)
  if (any(w != 0 & w != 1)) stop("w must be a 0/1 indicator")
  use <- w != 0
  Xu <- X[use, , drop = FALSE]; yu <- y[use]
  if (sum(yu) == 0 || sum(yu) == length(yu))
    stop(sprintf("the outcome has no %s among the %d subjects with mediator and outcome observed",
                 if (sum(yu) == 0) "events" else "non-events", length(yu)))
  qx <- qr(Xu, tol = 1e-7)
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  Xk <- Xu[, keep, drop = FALSE]
  b <- numeric(ncol(Xk)); conv <- FALSE
  for (it in seq_len(maxit)) {
    p <- plogis(drop(Xk %*% b)); W <- p * (1 - p)
    Iinv <- tryCatch(solve(crossprod(Xk * W, Xk)), error = function(e)
      stop("logistic fit failed (separation?)", if (!firth) "; use firth = TRUE", call. = FALSE))
    U <- if (firth) {
      h <- W * rowSums((Xk %*% Iinv) * Xk)
      crossprod(Xk, yu - p + h * (0.5 - p))
    } else crossprod(Xk, yu - p)
    step <- drop(Iinv %*% U)
    mx <- max(abs(step)); if (mx > 2) step <- step * 2 / mx
    b <- b + step
    if (mx < tol) { conv <- TRUE; break }
  }
  if (!conv) warning(sprintf("logistic fit did not converge in %d iterations", maxit), call. = FALSE)
  p <- plogis(drop(Xk %*% b)); W <- p * (1 - p)
  J <- crossprod(Xk * W, Xk) / n
  coef <- numeric(ncol(X)); coef[keep] <- b; names(coef) <- colnames(X)
  IF <- matrix(0, n, ncol(X)); colnames(IF) <- colnames(X)
  IF[use, keep] <- (Xk * (yu - p)) %*% solve(J)
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

.numjac <- function(fun, par) {
  h <- 1e-6 * (1 + abs(par))
  sapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, h[j])
    (fun(par + e) - fun(par - e)) / (2 * h[j])
  })
}

.fit_seq_binary <- function(Y, M, S, Z, Zt, binary_method = c("proj", "closed"),
                            firth = TRUE, gh_nodes = 15) {
  binary_method <- match.arg(binary_method)
  S <- as.matrix(S)
  n <- nrow(Z); p <- ncol(Z); k <- ncol(S)
  if (is.null(colnames(S)) || any(colnames(S) == "")) colnames(S) <- paste0("S", seq_len(k))
  R1 <- as.numeric(!is.na(M)); R2 <- R1 * !is.na(Y)
  Mz <- ifelse(R1 == 1, M, 0); Yz <- ifelse(R2 == 1, Y, 0)

  f1 <- lapply(seq_len(k), function(j) .wls_if(Z, S[, j]))
  kappa <- matrix(sapply(f1, `[[`, "coef"), p, k)
  E_S <- matrix(sapply(f1, `[[`, "resid"), n, k)
  Sigma_S <- crossprod(E_S) / n
  f2 <- .wls_if(cbind(S, Z), M, w = R1)
  sigma2_M <- sum(f2$resid^2) / sum(R1)
  f3 <- .logit_if(cbind(M = Mz, S, Z), Yz, w = R2, firth = firth)
  key <- c(colnames(S), colnames(Z)[2], "M")
  bad <- intersect(c(unlist(lapply(f1, `[[`, "aliased")), f2$aliased, f3$aliased), key)
  if (length(bad))
    stop("not identified (collinear within the observed subset): ", paste(bad, collapse = ", "))
  IF_s2 <- R1 * (f2$resid^2 - sigma2_M) * n / sum(R1)
  ia <- seq_len(p); ib <- p + seq_len(p + 1)

  if (binary_method == "closed") {
    lt <- lower.tri(Sigma_S, diag = TRUE)
    IF_Sig <- do.call(cbind, lapply(which(lt), function(ix) {
      j <- row(Sigma_S)[ix]; l <- col(Sigma_S)[ix]
      E_S[, j] * E_S[, l] - Sigma_S[j, l]
    }))
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
      vM <- drop(crossprod(delta, Sig %*% delta) + s2)
      nu <- drop(Sig %*% delta) / vM
      V <- Sig - tcrossprod(nu) * vM
      bS <- sum(omS * nu)
      bZ <- omZ + drop(kap %*% omS) - alpha * bS
      bZ[1] <- bZ[1] + 0.5 * drop(crossprod(omS, V %*% omS))
      c(alpha, bZ, omM + bS)
    }
    IF_par <- cbind(do.call(cbind, lapply(f1, `[[`, "IF")), IF_Sig, f2$IF, IF_s2, f3$IF)
    est <- map(par)
    IF <- IF_par %*% t(.numjac(map, par))
    alpha <- est[ia]; theta <- est[ib]
    IF_a <- IF[, ia, drop = FALSE]; IF_b <- IF[, ib, drop = FALSE]
    aliased <- f3$aliased
  } else {
    para <- c(kappa, f2$coef)
    mapa <- function(par) {
      kap <- matrix(par[seq_len(p * k)], p, k); dg <- par[-seq_len(p * k)]
      dg[-(1:k)] + drop(kap %*% dg[1:k])
    }
    alpha <- mapa(para)
    IF_a <- cbind(do.call(cbind, lapply(f1, `[[`, "IF")), f2$IF) %*% t(.numjac(mapa, para))
    gh <- .gh_rule(gh_nodes)
    o <- which(R1 == 1); u <- which(R1 == 0); n0 <- length(u)
    id <- c(o, rep(u, each = gh_nodes))
    wq <- c(rep(1, length(o)), rep(gh$w, times = n0))
    zq <- c(rep(0, length(o)), rep(gh$z, times = n0))
    isu <- c(rep(FALSE, length(o)), rep(TRUE, n0 * gh_nodes))
    yobs <- c(R2[o] == 1, rep(FALSE, n0 * gh_nodes))
    Zs <- Z[id, , drop = FALSE]; SZs <- cbind(S, Z)[id, , drop = FALSE]
    SZu <- SZs[isu, , drop = FALSE]; Mbase <- Mz[id]; Ybase <- Yz[id]
    q2 <- k + p; q3 <- 1 + k + p
    xi <- c(f2$coef, sigma2_M, f3$coef)
    stack <- function(xi) {
      dg <- xi[seq_len(q2)]; s2 <- xi[q2 + 1]; om <- xi[q2 + 1 + seq_len(q3)]
      m <- Mbase
      m[isu] <- drop(SZu %*% dg) + sqrt(s2) * zq[isu]
      r <- plogis(om[1] * m + drop(SZs %*% om[-1]))
      r[yobs] <- Ybase[yobs]
      list(x = cbind(Zs, M = m), r = r)
    }
    st <- stack(xi)
    th <- numeric(p + 1); conv <- FALSE
    for (it in 1:100) {
      pr <- plogis(drop(st$x %*% th)); W <- wq * pr * (1 - pr)
      step <- drop(solve(crossprod(st$x * W, st$x), crossprod(st$x, wq * (st$r - pr))))
      mx <- max(abs(step)); if (mx > 2) step <- step * 2 / mx
      th <- th + step
      if (mx < 1e-10) { conv <- TRUE; break }
    }
    if (!conv) warning("the final logistic regression did not converge", call. = FALSE)
    theta <- th
    pr <- plogis(drop(st$x %*% theta))
    ps <- rowsum(st$x * (wq * (st$r - pr)), id, reorder = FALSE)
    psi <- matrix(0, n, p + 1); psi[as.integer(rownames(ps)), ] <- ps
    A <- crossprod(st$x * (wq * pr * (1 - pr)), st$x) / n
    Psi <- function(xi) {
      s <- stack(xi)
      colSums(s$x * (wq * (s$r - plogis(drop(s$x %*% theta))))) / n
    }
    IF_xi <- cbind(f2$IF, IF_s2, f3$IF)
    IF_b <- (psi + IF_xi %*% t(.numjac(Psi, xi))) %*% t(solve(A))
    aliased <- character()
  }
  names(alpha) <- colnames(Z); colnames(IF_a) <- colnames(Z)
  names(theta) <- colnames(Zt); colnames(IF_b) <- colnames(Zt)
  rz <- .wls_if(Z, M, w = R1)$resid[R1 == 1]
  list(coef = theta, vcov = .vcov_if(IF_b), IF = IF_b,
       method = sprintf("SynMediation-logistic (%s%s)", binary_method, if (firth) ", Firth" else ""),
       mediator = list(coef = alpha, IF = IF_a, vcov = .vcov_if(IF_a), method = "SynSurr"),
       steps = list(omega = f3$coef, sigma2_M = sigma2_M, beta_M_direct = f3$coef[["M"]],
                    beta_M_surrogate = theta[[p + 1]] - f3$coef[["M"]],
                    partial_r2 = 1 - sigma2_M / mean(rz^2),
                    events = f3$events, step3_aliased = f3$aliased),
       aliased = aliased, n_obs_both = sum(R2))
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
