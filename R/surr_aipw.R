# Surrogate-augmented AIPW for the mediator-outcome path when both the mediator
# and the outcome are partially missing.
#
# The missingness is treated as monotone: level 0 observes W = (1, A, X, S_M,
# S_Y); level 1 adds M (R1 = 1); level 2 adds Y (R2 = R1 * R_Y = 1). Subjects
# with Y but no M are treated as level 0. The estimating function is the
# monotone-coarsening AIPW
#   R2/pi2 U(beta) + (1 - R1/pi1) q1(W) + (R1/pi1 - R2/pi2) q2(W, M),
# with U = Zt (Y - Zt' beta), pi1 = P(R1 = 1 | W), pi2 = pi1 P(R2 = 1 | R1 = 1, W, M),
# q1 = E[U | W] and q2 = E[U | W, M]. q2 predicts the missing outcome from its
# surrogate given the regressors, as in SynSurr; q1 integrates over the missing
# mediator. Working models: Y | W, M linear (fitted on R2 = 1) and
# M | W ~ N(W'lambda, sigma^2) (fitted on R1 = 1), so q1 and q2 are linear in
# beta and the estimator has a closed form.
.fit_surr_aipw <- function(Y, M, C, W, Zt, ps_trim = NULL) {
  n <- length(Y)
  pc <- ncol(C)
  K <- ncol(Zt)
  R1 <- as.numeric(!is.na(M))
  R2 <- R1 * !is.na(Y)
  Mz <- ifelse(R1 == 1, M, 0)
  Yz <- ifelse(R2 == 1, Y, 0)
  WM <- cbind(W, M = Mz)

  # Propensity of observing M, and of observing Y among those with M.
  ps1 <- .fit_propensity(W, R1, ps_trim)
  pi1 <- ps1$omega
  lam2 <- rep(1, n)
  score <- ps1$score
  if (any(R1 == 1 & R2 == 0)) {
    use <- R1 == 1
    ps2 <- .fit_propensity(WM[use, , drop = FALSE], R2[use], ps_trim)
    lam2[use] <- ps2$omega
    s2 <- matrix(0, n, ncol(WM))
    s2[use, ] <- ps2$score
    score <- cbind(score, s2)
  }
  pi2 <- pi1 * lam2

  # Working models.
  wm_M <- .wls_if(W, M, w = R1)
  m1 <- drop(W %*% wm_M$coef)
  m2 <- m1^2 + sum(wm_M$resid[R1 == 1]^2) / (sum(R1) - ncol(W))
  wm_Y <- .wls_if(WM, Y, w = R2)
  th_W <- wm_Y$coef[seq_len(ncol(W))]
  th_M <- wm_Y$coef[[ncol(WM)]]
  eyw <- drop(W %*% th_W)
  mu2 <- eyw + th_M * Mz                      # E[Y | W, M]
  EZY1 <- cbind(C * (eyw + th_M * m1), m1 * eyw + m2 * th_M)   # E[Zt Y | W]

  sum_EZZ <- function(c) {                    # sum_i c_i E[Zt_i Zt_i' | W_i]
    cm <- crossprod(C, c * m1)
    rbind(cbind(crossprod(C * c, C), cm), c(cm, sum(c * m2)))
  }
  a2 <- R2 / pi2
  a1 <- 1 - R1 / pi1
  b <- R1 / pi1 - R2 / pi2

  A_mat <- crossprod(Zt * (a2 + b), Zt) + sum_EZZ(a1)
  b_vec <- drop(crossprod(Zt * a2, Yz) + crossprod(EZY1, a1) + crossprod(Zt * b, mu2))
  beta <- drop(solve(A_mat, b_vec))
  names(beta) <- colnames(Zt)

  lc <- drop(C %*% beta[1:pc])
  q1 <- EZY1 - cbind(C * (lc + m1 * beta[K]), m1 * lc + m2 * beta[K])
  psi <- Zt * (a2 * (Yz - drop(Zt %*% beta)) + b * (mu2 - drop(Zt %*% beta))) + a1 * q1

  J <- A_mat / n
  IF <- .residualize(psi, score) %*% t(solve(J))
  colnames(IF) <- colnames(Zt)
  list(coef = beta, vcov = .vcov_if(IF), IF = IF, method = "SurrAIPW",
       pi1 = pi1, pi2 = pi2, working_model_M = wm_M$coef,
       working_model_Y = wm_Y$coef, n_obs_both = sum(R2))
}
