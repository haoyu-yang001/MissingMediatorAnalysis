# Base-R checks run by R CMD check.
library(MissingMediatorAnalysis)

set.seed(1)
d <- sim_mma_data(n = 3000, rho = 0.8, miss_rate = 0.5, mechanism = "MAR")
X <- cbind(1, d$A, d$X)

# 1. Without missingness both paths reduce to ordinary least squares.
d_full <- transform(d, M = M_full)
f0 <- mma(d_full, "A", "M", "Y", "S", "X")
stopifnot(abs(f0$effects["alpha", "estimate"] - coef(lm(M_full ~ A + X, d))[["A"]]) < 1e-10)
stopifnot(abs(f0$effects["beta", "estimate"] - coef(lm(Y ~ A + X + M_full, d))[["M_full"]]) < 1e-10)

# 2. Complete-case standard errors equal the HC0 sandwich.
cc <- mma(d, "A", "M", "Y", "S", "X", mediator_method = "CC", outcome_method = "CC")
dc <- d[!is.na(d$M), ]
l <- lm(M ~ A + X, dc)
Xc <- model.matrix(l)
B <- solve(crossprod(Xc))
hc0 <- B %*% crossprod(Xc * resid(l)) %*% B
stopifnot(abs(cc$effects["alpha", "se"] - sqrt(hc0["A", "A"])) < 1e-10)

# 3. SynSurr equals the two-step closed form gamma + delta * alpha_S.
ss <- fit_synsurr(d$M, d$S, X)
a_S <- coef(lm(S ~ A + X, d))
b <- coef(lm(M ~ S + A + X, d))
stopifnot(max(abs(ss$coef - (b[c("(Intercept)", "A", "X")] + b[["S"]] * a_S))) < 1e-10)

# 4. AIPW_new converges and its point estimate is close to AIPW when the
#    working model is correct.
W <- cbind(1, d$Y, d$A, d$X, d$S)
o1 <- fit_outcome_model(d$Y, d$M, X, W, "AIPW")
o2 <- fit_outcome_model(d$Y, d$M, X, W, "AIPWnew")
stopifnot(o2$n_iter < 100, max(abs(o1$coef - o2$coef)) < 0.01)

# 5. The indirect-effect variance from stacked influence functions matches
#    the delta method with the estimated cross-path covariance.
f <- mma(d, "A", "M", "Y", "S", "X")
V <- f$vcov
a <- f$effects["alpha", "estimate"]; bM <- f$effects["beta", "estimate"]
v_delta <- bM^2 * V["alpha", "alpha"] + a^2 * V["beta", "beta"] + 2 * a * bM * V["alpha", "beta"]
stopifnot(abs(v_delta - V["NIE", "NIE"]) < 1e-12)

# 6. The combined estimator is no less efficient than either component, and its
#    variance equals the closed-form minimum-variance combination.
oc <- fit_outcome_model(d$Y, d$M, X, W, "Combined")
v <- function(f) diag(f$vcov)
stopifnot(all(v(oc) <= v(oc$components$AIPWnew) + 1e-12),
          all(v(oc) <= v(oc$components$CC) + 1e-12))
I1 <- oc$components$AIPWnew$IF[, 3]; I2 <- oc$components$CC$IF[, 3]
S <- crossprod(cbind(I1, I2)) / nrow(d)^2
v_opt <- (S[1, 1] * S[2, 2] - S[1, 2]^2) / (S[1, 1] + S[2, 2] - 2 * S[1, 2])
stopifnot(v(oc)[3] <= v_opt + 1e-12)

# 7. Monte Carlo NIE interval: reproducible under set.seed, close to the delta
#    interval when both paths are far from zero, and skewed otherwise.
set.seed(2); m1 <- mma(d, "A", "M", "Y", "S", "X")
set.seed(2); m2 <- mma(d, "A", "M", "Y", "S", "X")
stopifnot(identical(m1$nie_ci, m2$nie_ci), nrow(m1$nie_ci) == 2)
stopifnot(nrow(mma(d, "A", "M", "Y", "S", "X", mc_draws = 0)$nie_ci) == 1)
ci <- MissingMediatorAnalysis:::.mc_product_ci(c(5, 5), diag(2) * 0.01, 2e5)
stopifnot(abs(ci[1] - (25 - 1.96 * sqrt(0.5))) < 0.02,
          abs(ci[2] - (25 + 1.96 * sqrt(0.5))) < 0.02)
ci0 <- MissingMediatorAnalysis:::.mc_product_ci(c(0.1, 0.1), diag(2) * 0.05^2, 2e5)
stopifnot(ci0[2] - 0.01 > 0.01 - ci0[1])

# 8. Missing outcome: with complete data sequential SynSurr reproduces least
#    squares; with a missing outcome mma() defaults to it and stays near the
#    truth in a large sample with missingness that depends on both surrogates.
set.seed(3)
dy <- sim_mma_data(n = 3000, rho = 0.8, rho_Y = 0.8, miss_rate = 0, miss_rate_Y = 0)
fs <- mma(transform(dy, M = M_full, Y = Y_full), "A", "M", "Y", "S", "X",
          outcome_surrogate = "S_Y", outcome_method = "SeqSynSurr", mc_draws = 0)
stopifnot(abs(fs$effects["beta", "estimate"] - coef(lm(Y_full ~ A + X + M_full, dy))[["M_full"]]) < 1e-10,
          abs(fs$effects["alpha", "estimate"] - coef(lm(M_full ~ A + X, dy))[["A"]]) < 1e-10)
dy <- sim_mma_data(n = 1e5, rho = 0.8, miss_rate = 0.8, mechanism = "MAR",
                   rho_Y = 0.8, miss_rate_Y = 0.5, mechanism_Y = "MAR")
fy <- mma(dy, "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y", mc_draws = 0)
stopifnot(fy$methods[["outcome"]] == "SeqSynSurr", fy$n_obs_Y < fy$n,
          all(abs(fy$effects[c("alpha", "beta", "NDE"), "estimate"] - 0.1) <
                4 * fy$effects[c("alpha", "beta", "NDE"), "se"]))

# 9. With a complete outcome and no outcome surrogate the default is sequential
#    SynSurr with the mediator surrogate only.
f9 <- mma(d, "A", "M", "Y", "S", "X", mc_draws = 0)
stopifnot(f9$methods[["outcome"]] == "SeqSynSurr",
          abs(f9$effects["alpha", "estimate"] - fit_synsurr(d$M, d$S, X)$coef[[2]]) < 1e-10)

# 10. Robustness for large-scale use: an indicator covariate that is zero for
#     every subject with the mediator is reported as NA without affecting the
#     path estimates; factor covariates are expanded; rows with a missing
#     exposure are dropped; the exact Monte Carlo p-value agrees with draws.
set.seed(4)
d10 <- sim_mma_data(n = 4000, rho = 0.8, miss_rate = 0.6, rho_Y = 0.8, miss_rate_Y = 0.3)
d10$ctr <- factor(sample(c("a", "b", "c"), 4000, TRUE))
d10$ctr[is.na(d10$M)][1:50] <- "c"; d10$ctr[!is.na(d10$M) & d10$ctr == "c"] <- "a"
d10$Xs <- d10$X
f10 <- mma(d10, "A", "M", "Y", "S", c("X", "ctr"), outcome_surrogate = "S_Y", mc_draws = 0)
f10b <- mma(d10[d10$ctr != "c", ], "A", "M", "Y", "S", c("X", "ctr"), outcome_surrogate = "S_Y", mc_draws = 0)
stopifnot("ctrc" %in% f10$aliased, is.na(f10$mediator_model["ctrc", "estimate"]),
          abs(f10$effects["alpha", "estimate"] - f10b$effects["alpha", "estimate"]) < 0.02)
d10$A[1:10] <- NA
f10c <- suppressMessages(mma(d10, "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y", mc_draws = 0))
stopifnot(f10c$n == 3990)
V <- matrix(c(0.01, 0.003, 0.003, 0.02), 2)
z <- matrix(rnorm(4e6), ncol = 2) %*% chol(V)
p_draw <- mean((0.15 + z[, 1]) * (0.1 + z[, 2]) < 0)
stopifnot(abs(MissingMediatorAnalysis:::.product_p(c(0.15, 0.1), V) - 2 * min(p_draw, 1 - p_draw)) < 0.003)

# 11. mma_scan returns one row per pair and records failures; diagnostics run.
G <- cbind(A = d10$A, junk = c(rnorm(8), rep(NA, 3992)))
sc <- mma_scan(d10, G, mediators = "M", surrogates = "S", outcome = "Y",
               covariates = "X", outcome_surrogate = "S_Y", mc_draws = 0)
stopifnot(nrow(sc) == 2, sc$status[2] == "error", sc$status[1] == "ok")
dg <- mma_diagnostics(f10, d10[!is.na(d10$A), ])
stopifnot(all(c("linearity", "variance", "missingness", "hausman") %in% names(dg)))

# 12. Auxiliary variables: an extra fully observed column in the surrogate set
#     (missingness depends on it) keeps sequential SynSurr near the truth, and
#     the complete-data fit still equals least squares.
set.seed(5)
d12 <- sim_mma_data(n = 2000, rho = 0.8, rho_Y = 0.8, miss_rate = 0, miss_rate_Y = 0)
d12$aux <- rnorm(2000)
f12 <- mma(transform(d12, M = M_full, Y = Y_full), "A", "M", "Y", c("S", "aux"), "X",
           outcome_surrogate = "S_Y", mc_draws = 0)
stopifnot(abs(f12$effects["beta", "estimate"] - coef(lm(Y_full ~ A + X + M_full, d12))[["M_full"]]) < 1e-10)
err <- tryCatch(mma(d12, "A", "M", "Y", c("S", "aux"), "X", outcome_method = "Combined"),
                error = function(e) "error")
stopifnot(identical(err, "error"))
# diagnostics with an auxiliary column (complete-case refit uses S only) and the
# per-column table
f12m <- mma(d12, "A", "M", "Y", c("S", "aux"), "X", outcome_surrogate = "S_Y", mc_draws = 0)
dg12 <- mma_diagnostics(f12m, d12)
stopifnot(identical(dg12$surrogate_columns$column, c("S", "aux", "S_Y")),
          is.na(dg12$surrogate_columns$beta_without[1]),
          all(is.finite(dg12$surrogate_columns$se_ratio_with_vs_without[-1])))

# 13. Nonlinear step 2 (basis): equals least squares with complete data, runs
#     with missing data, and diagnostics report the heteroskedasticity and
#     nonlinearity checks.
d13 <- transform(d12, M = M_full, Y = Y_full)
ols <- coef(lm(Y_full ~ A + X + M_full, d12))[["M_full"]]
for (b in c("poly2", "ns3")) {
  f13 <- mma(d13, "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y", mc_draws = 0, basis = b)
  stopifnot(abs(f13$effects["beta", "estimate"] - ols) < 1e-10,
            abs(f13$effects["alpha", "estimate"] - coef(lm(M_full ~ A + X, d12))[["A"]]) < 1e-10)
}
f13 <- mma(d10[!is.na(d10$A), ], "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y",
           mc_draws = 0, basis = "poly2")
stopifnot(f13$methods[["outcome"]] == "SeqSynSurr(poly2)", all(is.finite(f13$effects$se[1:5])))
stopifnot(all(c("heteroskedasticity", "nonlinearity") %in% names(dg)))
f13b <- mma(d13, "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y", mc_draws = 0,
            basis = "ns3", basis_outcome = TRUE)
stopifnot(abs(f13b$effects["beta", "estimate"] - ols) < 1e-10)

# 14. Sequential SynSurr without an outcome surrogate (missing outcome).
f14 <- mma(d10[!is.na(d10$A), ], "A", "M", "Y", "S", "X", mc_draws = 0)
stopifnot(f14$methods[["outcome"]] == "SeqSynSurr", is.finite(f14$effects["beta", "se"]))

# 15. mma_scan with a seed is reproducible.
s1 <- mma_scan(d10, G[, 1, drop = FALSE], "M", "S", "Y", "X", "S_Y", mc_draws = 1000, seed = 1)
s2 <- mma_scan(d10, G[, 1, drop = FALSE], "M", "S", "Y", "X", "S_Y", mc_draws = 1000, seed = 1)
stopifnot(identical(s1$mc_lo, s2$mc_lo))

cat("all tests passed\n")

# Binary outcome: target-preserving sequential AIPW / IPW (logistic outcome model).
set.seed(7)
nb <- 6000
gb <- data.frame(A = rnorm(nb), X = rnorm(nb))
gb$M_full <- 0.3 * gb$A + 0.2 * gb$X + rnorm(nb)
gb$S <- gb$M_full + rnorm(nb)
gb$S_Y <- rnorm(nb)
gb$Y_full <- rbinom(nb, 1, plogis(-2 + 0.5 * gb$M_full + 0.2 * gb$A + 0.2 * gb$X + 0.8 * gb$S_Y))
miss_M <- function(eta) ifelse(rbinom(nb, 1, plogis(eta)) == 1, gb$M_full, NA)
# (a) No missing data: AIPW (any firth) and IPW with firth = FALSE reproduce lm / glm.
dfull <- transform(gb, M = M_full, Y = Y_full)
g0 <- glm(Y ~ A + X + M, binomial, dfull)
for (om in c("SeqSynSurr", "IPW")) {
  f <- mma(dfull, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_method = om,
           firth = om != "IPW", mc_draws = 0)
  stopifnot(max(abs(f$fits$outcome$coef - coef(g0))) < 1e-7,
            abs(f$effects["alpha", "estimate"] - coef(lm(M_full ~ A + X, gb))[["A"]]) < 1e-10,
            f$missing_data$order == "none")
}
# (b) Complete-case ML logistic equals glm on complete cases.
db <- transform(gb, M = miss_M(-0.5 + 0.5 * S + 0.5 * Y_full), Y = Y_full)
fcc <- mma(db, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_method = "CC",
           mediator_method = "CC", firth = FALSE, mc_draws = 0)
gcc <- glm(Y ~ A + X + M, binomial, db[!is.na(db$M), ])
stopifnot(max(abs(fcc$fits$outcome$coef - coef(gcc))) < 1e-7)
# (c) Y complete, M missing at random given (S, A, X, Y): IPW with firth = FALSE equals a
#     hand-coded inverse-probability-weighted fit, with the influence function
#     residualised on the missingness-model score; the default (sequential AIPW) uses
#     the ordering Y -> M with the outcome in the missingness model.
fi <- mma(db, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_method = "IPW",
          firth = FALSE, mc_draws = 0)
R <- !is.na(db$M)
V <- cbind(1, db$S, db$A, db$X, db$Y)
ps <- glm(R ~ V - 1, family = binomial)
w <- R / fitted(ps)
dw <- db[R, ]; dw$w <- w[R]
gw <- suppressWarnings(glm(Y ~ A + X + M, quasibinomial, dw, weights = w))
lw <- lm(M ~ A + X, dw, weights = w)
xm <- model.matrix(gw); mu <- fitted(gw)
J <- crossprod(xm * (dw$w * mu * (1 - mu)), xm) / nb
IF <- matrix(0, nb, 4); IF[R, ] <- (xm * (dw$w * (dw$Y - mu))) %*% solve(J)
sc <- (R - fitted(ps)) * V
IFr <- IF - sc %*% solve(crossprod(sc), crossprod(sc, IF))
stopifnot(max(abs(fi$fits$outcome$coef - coef(gw))) < 1e-7,
          abs(fi$effects["alpha", "estimate"] - coef(lw)[["A"]]) < 1e-10,
          abs(fi$effects["beta", "se"] - sqrt(sum(IFr[, 4]^2)) / nb) < 1e-8,
          fi$missing_data$order == "YM")
fa <- mma(db, "A", "M", "Y", "S", "X", outcome_family = "binomial", mc_draws = 0)
d_a <- fa$binary_diagnostics
stopifnot(fa$missing_data$order == "YM", grepl("AIPW", fa$methods[["outcome"]]),
          isTRUE(d_a$converged), d_a$eq_norm < 1e-10, d_a$n_MY == sum(R),
          abs(d_a$max_weight - max(w)) < 1e-8,
          abs(fa$effects["beta", "estimate"] - 0.5) < 4 * fa$effects["beta", "se"],
          names(fa$fits$outcome$coef)[4] == "M", nrow(fa$effects_or) == 4)
# Stacked influence functions give the delta-method NIE variance.
Vb <- fa$vcov; a <- fa$effects["alpha", "estimate"]; bM <- fa$effects["beta", "estimate"]
stopifnot(abs(bM^2 * Vb["alpha", "alpha"] + a^2 * Vb["beta", "beta"] + 2 * a * bM * Vb["alpha", "beta"] -
              Vb["NIE", "NIE"]) < 1e-12)
# (d) Missing-data ordering: monotone -> M -> Y; outcome without mediator -> Y -> M;
#     missing_order overrides; ps_trim bounds the weights.
dm <- transform(gb, M = miss_M(0.5 * S))
dm$Y <- ifelse(is.na(dm$M) | rbinom(nb, 1, 0.8) == 1, ifelse(is.na(dm$M), NA, dm$Y_full), NA)
f_my <- mma(dm, "A", "M", "Y", "S", "X", outcome_family = "binomial", mc_draws = 0)
dym <- transform(gb, M = miss_M(0.5 * S), Y = ifelse(rbinom(nb, 1, 0.8) == 1, Y_full, NA))
f_ym <- mma(dym, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_surrogate = "S_Y",
            mc_draws = 0)
f_ym2 <- mma(dym, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_surrogate = "S_Y",
             missing_order = "MY", mc_draws = 0)
f_tr <- mma(dym, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_surrogate = "S_Y",
            ps_trim = 0.5, mc_draws = 0)
stopifnot(f_my$missing_data$order == "MY", f_ym$missing_data$order == "YM",
          f_ym2$missing_data$order == "MY", f_ym2$binary_diagnostics$n_MY == f_ym$binary_diagnostics$n_MY,
          f_tr$binary_diagnostics$max_weight <= 1 / 0.5^2 + 1e-12,
          all(is.finite(f_my$effects$se[1:5])), all(is.finite(f_ym$effects$se[1:5])))
# (e) mma_scan with a binary outcome; keep_IF = FALSE; invalid inputs are rejected.
scb <- mma_scan(dym, c("A", "X"), "M", "S", "Y", outcome_surrogate = "S_Y",
                outcome_family = "binomial", mc_draws = 0)
stopifnot(nrow(scb) == 2, all(scb$status != "error"), all(scb$missing_order == "YM"),
          all(is.finite(scb$max_if_share_beta)))
fk <- mma(dym, "A", "M", "Y", "S", "X", outcome_family = "binomial", outcome_surrogate = "S_Y",
          keep_IF = FALSE, mc_draws = 0)
stopifnot(is.null(fk$fits$outcome$IF), is.null(fk$fits$mediator$IF))
invisible(capture.output(print(f_ym), summary(f_my)))
bad <- function(...) inherits(try(mma(dym, "A", "M", ..., outcome_family = "binomial"),
                                  silent = TRUE), "try-error")
stopifnot(bad("M_full", "S", "X"), bad("Y", "S", "X", outcome_method = "Combined"),
          bad("Y", "S", "X", mediator_method = "CC"))
cat("binary tests passed\n")
