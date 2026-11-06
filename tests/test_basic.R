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

# Binary outcome (logistic outcome model).
set.seed(7)
db <- sim_mma_data(n = 4000, miss_rate = 0.6, rho = 0.8, rho_Y = 0.6, miss_rate_Y = 0.4)
db$Yb <- ifelse(is.na(db$Y), NA, as.numeric(db$Y > quantile(db$Y, 0.85, na.rm = TRUE)))
# (a) With no missing data and firth = FALSE, proj reproduces glm.
dfull <- transform(db, M = M_full)
dfull$Yb <- as.numeric(dfull$Y_full > quantile(dfull$Y_full, 0.85))
fb0 <- mma(dfull, "A", "M", "Yb", "S", "X", outcome_family = "binomial", firth = FALSE, mc_draws = 0)
g0 <- glm(Yb ~ A + X + M, binomial, dfull)
stopifnot(abs(fb0$effects["beta", "estimate"] - coef(g0)[["M"]]) < 1e-7)
stopifnot(abs(fb0$effects["NDE", "estimate"] - coef(g0)[["A"]]) < 1e-7)
# (b) Complete-case ML logistic equals glm on complete cases.
fcc <- mma(db, "A", "M", "Yb", "S", "X", outcome_family = "binomial", outcome_method = "CC",
           mediator_method = "CC", firth = FALSE, mc_draws = 0)
gcc <- glm(Yb ~ A + X + M, binomial, db[!is.na(db$M) & !is.na(db$Yb), ])
stopifnot(abs(fcc$effects["beta", "estimate"] - coef(gcc)[["M"]]) < 1e-7)
# (c) Proposed fit runs, proj and closed are close for this outcome, OR table present.
fb <- mma(db, "A", "M", "Yb", "S", "X", outcome_surrogate = "S_Y", outcome_family = "binomial", mc_draws = 0)
fbc <- mma(db, "A", "M", "Yb", "S", "X", outcome_surrogate = "S_Y", outcome_family = "binomial",
           binary_method = "closed", mc_draws = 0)
stopifnot(is.finite(fb$effects["NIE", "se"]), nrow(fb$effects_or) == 4,
          abs(fb$effects["beta", "estimate"] - fbc$effects["beta", "estimate"]) < 3 * fb$effects["beta", "se"])
# (d) Non-binary outcome is rejected.
stopifnot(inherits(try(mma(db, "A", "M", "Y", "S", "X", outcome_surrogate = "S_Y",
                           outcome_family = "binomial"), silent = TRUE), "try-error"))
cat("binary tests passed\n")
