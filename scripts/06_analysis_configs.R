## Tables 6, 7, 8, 12: mean/variance ratios, skewness, lattice Kolmogorov distances, prediction-interval coverage
source("scripts/00_setup.R")
Z <- fromJSON(file.path(OUT, "zeta_configs.json"))
skew <- function(x) mean((x - mean(x))^3)/mean((x - mean(x))^2)^1.5   # moment (g1) skewness
exkurt <- function(x) mean((x - mean(x))^4)/mean((x - mean(x))^2)^2 - 3
lattice_dk <- function(x, G) {      # sup over integers of |F_hat(k) - G(k)| (Chun & Jhong 2026b)
  x <- sort(x); ks <- (x[1] - 1):x[length(x)]
  F <- findInterval(ks, x)/length(x)
  max(abs(F - G(ks)), G(x[1] - 1), 1 - G(x[length(x)]))
}
model <- function(nm) {
  cf <- CFG[[nm]]; Lb <- sum(cf$lens*cf$probs); varL <- sum(cf$probs*cf$lens^2) - Lb^2; r <- cf$r
  if (isTRUE(cf$torus)) return(c(q = pi*r*Lb^2 + 8*pi*r^2*Lb + 32/3*pi*r^3, zeta = (pi*r*(Lb + 4*r))^2*varL))
  c <- cf$c %||% 1; Es <- if (c == 1) pi/4 else Z[[paste0("Esin_", c)]]$mean
  c(q = q_exact(r, cf$lens, cf$probs, RHO, cf$box %||% c(1,1,1), Es), zeta = Z[[nm]]$zeta)
}
out <- NULL
for (nm in names(CFG)) {
  fn <- file.path(OUT, sprintf("counts_%s.csv", nm)); if (!file.exists(fn)) next
  A <- as.matrix(read.csv(fn)); m <- model(nm); q <- m[["q"]]; z <- m[["zeta"]]
  for (g in seq_along(NG)) {
    n <- NG[g]; x <- A[, g]; R <- length(x); B <- n*(n-1)/2
    EN <- B*q; VN <- B*(2*(n-2)*z + q*(1-q)); sdN <- sqrt(VN); nstar <- if (z > 0) q/(2*z) else Inf
    row <- data.frame(config = nm, n = n, R = R, nstar = nstar, n_over_nstar = n/nstar,
      mean = mean(x), EN = EN, mean_ratio = mean(x)/EN, mean_ratio_se = sd(x)/sqrt(R)/EN,
      var = var(x), VN = VN, var_ratio = var(x)/VN, proj_frac = n*(n-1)*(n-2)*z/VN,
      skew = skew(x), exkurt = exkurt(x),
      dK_theory = lattice_dk(x, function(k) pnorm((k + .5 - EN)/sdN)),
      dK_sample = lattice_dk(x, function(k) pnorm((k + .5 - mean(x))/sd(x))))
    for (lev in c(.90, .95, .99)) row[[sprintf("cov_normal_%d", 100*lev)]] <- mean(abs(x - EN) <= qnorm(.5 + lev/2)*sdN + .5)
    if (VN > EN) {                       # moment-matched negative binomial (eq. 33)
      p <- EN/VN; size <- EN^2/(VN - EN)
      row$dK_NB <- lattice_dk(x, function(k) pnbinom(k, size, p))
      for (lev in c(.90, .95, .99)) row[[sprintf("cov_NB_%d", 100*lev)]] <-
        mean(x >= qnbinom((1 - lev)/2, size, p) & x <= qnbinom(.5 + lev/2, size, p))
    }
    out <- if (is.null(out)) row else merge(out, row, all = TRUE, sort = FALSE)
  }
}
write.csv(out, file.path(OUT, "config_analysis.csv"), row.names = FALSE)
print(out[, c("config","n","n_over_nstar","mean_ratio","var_ratio","skew","dK_theory","cov_normal_95")], digits = 3)
