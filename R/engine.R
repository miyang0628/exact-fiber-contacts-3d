## engine.R -- R port of engine.py: sampling, closed forms, and C++ kernels
## Usage: source("R/engine.R") from the project root (compiles R/engine.cpp once).
suppressPackageStartupMessages(library(Rcpp))
.here <- if (exists("ROOT")) ROOT else "."
sourceCpp(file.path(.here, "R", "engine.cpp"), cacheDir = file.path(.here, ".rcpp_cache"))

## directions: isotropic, vertical component scaled by c (oblate for c<1)
directions <- function(n, c = 1) {
  v <- matrix(rnorm(3 * n), n, 3); v[, 3] <- v[, 3] * c
  v / sqrt(rowSums(v^2))
}

## sample n fibers: midpoints uniform in S = C (+) cube(lmax/2 + rho), or on the unit torus
sample_fibers <- function(n, lens, probs, rho, box = c(1, 1, 1), c = 1, torus = FALSE) {
  L <- if (length(lens) == 1) rep(lens, n) else sample(lens, n, replace = TRUE, prob = probs)
  T <- directions(n, c)
  if (torus) X <- matrix(runif(3 * n), n, 3) else {
    m <- max(lens) / 2 + rho
    X <- sapply(1:3, function(k) runif(n, -m, box[k] + m))
  }
  list(X = X, T = T, L = L)
}

Svol <- function(lens, rho, box = c(1, 1, 1)) prod(box + 2 * (max(lens) / 2 + rho))

## exact contact probability (Theorem 4); Esin = E sin(gamma) of the direction law
q_exact <- function(r, lens, probs, rho, box = c(1, 1, 1), Esin = pi / 4) {
  Lb <- sum(probs * lens); S <- Svol(lens, rho, box)
  prod(box) * (4 * r * Lb^2 * Esin + 8 * pi * r^2 * Lb + 32 / 3 * pi * r^3) / S^2
}
q_lead <- function(r, lens, probs, rho, box = c(1, 1, 1), Esin = pi / 4) {
  Lb <- sum(probs * lens); prod(box) * 4 * r * Lb^2 * Esin / Svol(lens, rho, box)^2
}

## truncated length moments lambda1, lambda2 (box, isotropic; eq. 10-12)
lam_moments <- function(lens, probs, rho, box = c(1, 1, 1)) {
  M <- function(k) sum(probs * lens^k); a <- box; C <- prod(a); S <- Svol(lens, rho, box)
  s1 <- a[2] * a[3] + a[1] * a[3] + a[1] * a[2]; s2 <- sum(a)
  c(lambda1 = M(1) * C / S,
    lambda2 = (C * M(2) - s1 / 6 * M(3) + s2 / (9 * pi) * M(4) - M(5) / (40 * pi)) / S)
}

## tip moments E T, E T^2, E[lC T] (eq. 13; box version of the same expansion)
tip_moments <- function(lens, probs, rho, box = c(1, 1, 1)) {
  M <- function(k) sum(probs * lens^k); a <- box; C <- prod(a); S <- Svol(lens, rho, box)
  s1 <- a[2] * a[3] + a[1] * a[3] + a[1] * a[2]; s2 <- sum(a)
  ET <- 2 * C / S
  both <- (C - 0.5 * s1 * M(1) + 2 / (3 * pi) * s2 * M(2) - M(3) / (4 * pi)) / S
  one <- (C * M(1) - 0.25 * s1 * M(2) + 2 / (9 * pi) * s2 * M(3) - M(4) / (16 * pi)) / S
  c(ET = ET, ET2 = ET + 2 * both, ElT = 2 * one)
}

## second-order projection variance (Proposition 5) and leading order
zeta1_second <- function(r, lens, probs, rho, box = c(1, 1, 1)) {
  Lb <- sum(probs * lens); S <- Svol(lens, rho, box)
  l <- lam_moments(lens, probs, rho, box); t <- tip_moments(lens, probs, rho, box)
  VarL <- l[2] - l[1]^2; VarT <- t[2] - t[1]^2; Cv <- t[3] - l[1] * t[1]
  al <- pi * r * (Lb + 4 * r); be <- 2 * pi * r^2 * Lb
  unname((al^2 * VarL + 2 * al * be * Cv + be^2 * VarT) / S^2)
}
zeta1_lead <- function(r, lens, probs, rho, box = c(1, 1, 1)) {
  Lb <- sum(probs * lens); S <- Svol(lens, rho, box); l <- lam_moments(lens, probs, rho, box)
  unname((pi * r * Lb / S)^2 * (l[2] - l[1]^2))
}

## exact zeta_{1,r} by nested Monte Carlo (Appendix D)
zeta_mc <- function(r, lens, probs, rho, box = c(1, 1, 1), c = 1, K = 1e5, M = 200, seed = 1) {
  set.seed(seed)
  f <- sample_fibers(K, lens, probs, rho, box, c); S <- Svol(lens, rho, box)
  cum <- cumsum(probs); cum[length(cum)] <- 1
  est <- t(vapply(seq_len(K), function(i)
    h1_est(f$X[i, ], f$T[i, ], f$L[i], lens, cum, r, M, box, c), numeric(2)))
  h <- est[, 1] / S; vv <- est[, 2] / S^2
  zeta <- var(h) - mean(vv)
  se <- sqrt(var((h - mean(h))^2 - vv) / K)
  c(zeta = zeta, se = se, hmean = mean(h))
}

Esin_mc <- function(c, N = 4e6, seed = 7) {
  set.seed(seed); a <- directions(N, c); b <- directions(N, c)
  cr <- cbind(a[, 2] * b[, 3] - a[, 3] * b[, 2], a[, 3] * b[, 1] - a[, 1] * b[, 3], a[, 1] * b[, 2] - a[, 2] * b[, 1])
  s <- sqrt(rowSums(cr^2)); c(mean = mean(s), se = sd(s) / sqrt(N))
}

## cube-truncated length lC (vectorized)
lC <- function(X, T, L, box = c(1, 1, 1)) {
  lo <- -L / 2; hi <- L / 2
  for (k in 1:3) {
    s0 <- (0 - X[, k]) / T[, k]; s1 <- (box[k] - X[, k]) / T[, k]
    lo <- pmax(lo, pmin(s0, s1)); hi <- pmin(hi, pmax(s0, s1))
  }
  pmax(hi - lo, 0)
}

## experiment configurations (Table 6)
CFG <- list(
  base   = list(lens = c(.2, .4), probs = c(.5, .5), r = .01),
  r005   = list(lens = c(.2, .4), probs = c(.5, .5), r = .005),
  r002   = list(lens = c(.2, .4), probs = c(.5, .5), r = .002),
  r001   = list(lens = c(.2, .4), probs = c(.5, .5), r = .001),
  equal  = list(lens = .3, probs = 1, r = .01),
  p02    = list(lens = c(.2, .4), probs = c(.2, .8), r = .01),
  p08    = list(lens = c(.2, .4), probs = c(.8, .2), r = .01),
  three  = list(lens = c(.1, .3, .5), probs = rep(1 / 3, 3), r = .01),
  c075   = list(lens = c(.2, .4), probs = c(.5, .5), r = .01, c = .75),
  c050   = list(lens = c(.2, .4), probs = c(.5, .5), r = .01, c = .5),
  c025   = list(lens = c(.2, .4), probs = c(.5, .5), r = .01, c = .25),
  box    = list(lens = c(.2, .4), probs = c(.5, .5), r = .01, box = c(2, 1, .5)),
  torus_eq  = list(lens = .3, probs = 1, r = .01, torus = TRUE),
  torus_mix = list(lens = c(.2, .4), probs = c(.5, .5), r = .01, torus = TRUE)
)
NG <- c(100L, 200L, 400L, 800L, 1600L, 3200L)
RHO <- 0.01
`%||%` <- function(a, b) if (is.null(a)) b else a
