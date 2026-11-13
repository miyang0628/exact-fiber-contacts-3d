## Tables 9-10: radius estimators (exact, first-order, thin-fiber) and Wald / data-based / test-inversion intervals
source("scripts/00_setup.R")
source(file.path(ROOT, "R", "inference.R"))
NGe <- c(200L, 400L, 800L, 1600L, 3200L); R <- reps(2000, 60); z <- qnorm(.975)
rows <- NULL
for (r in c(.002, .005, .01)) {
  res <- par_reps(R, function() {
    f <- sample_fibers(max(NGe), LENS, PROBS, RHO)
    P <- pairs_box(f$X, f$T, f$L, r, c(1,1,1))
    t(sapply(NGe, function(n) {
      Pn <- P[P[, 2] < n, , drop = FALSE]; N <- nrow(Pn); B <- n*(n-1)/2; U <- N/B
      deg <- tabulate(c(Pn) + 1, nbins = n)
      s2 <- var(deg/(n-1)); zh <- max((s2 - U*(1-U)/(n-1))*(n-1)/(n-2), 0)
      VNd <- B*(2*(n-2)*zh + U*(1-U))
      rex <- inv_ex(U); r1 <- inv_first(U); rl <- inv_lead(U)
      se_mod <- sqrt(varN(rex, n))/(B*dqex(rex)); se_dat <- sqrt(VNd)/(B*dqex(rex))
      zold <- zeta1_lead(r1, LENS, PROBS, RHO)*(1+8*r1/LBAR)^2; q1 <- A1*r1*(1+8*r1/LBAR)
      se_old <- sqrt(B*(2*(n-2)*zold + q1*(1-q1)))/(B*A1*(1+16*r1/LBAR))
      ti <- ti_interval(N, n, z)
      c(r = r, n = n, N = N, rex = rex, rlead = rl, rfirst = r1, se_mod = se_mod, se_dat = se_dat,
        se_old = se_old, ti_lo = ti[1], ti_hi = ti[2])
    }))
  }, seed = round(1e6*r))
  rows <- rbind(rows, do.call(rbind, res)); cat("r =", r, "done\n")
}
write.csv(rows, file.path(OUT, "estimation_replicates.csv"), row.names = FALSE)
## summary (Tables 9-10)
D <- as.data.frame(rows); S <- NULL
for (r in unique(D$r)) for (n in NGe) {
  d <- D[D$r == r & D$n == n, ]; Rr <- nrow(d); p <- qex(r)
  e <- function(x) c(bias = mean(x - r), bias_se = sd(x - r)/sqrt(Rr), sd = sd(x), rmse = sqrt(mean((x - r)^2)))
  cov <- function(lo, hi) c(mean(lo <= r & r <= hi), mean(hi - lo))
  S <- rbind(S, data.frame(r = r, n = n, R = Rr, t(e(d$rex)), first_bias = mean(d$rfirst - r), first_bias_se = sd(d$rfirst)/sqrt(Rr),
    first_pseudo = inv_first(p) - r, lead_bias = mean(d$rlead - r), lead_pseudo = inv_lead(p) - r,
    se_mod = mean(d$se_mod), se_dat = mean(d$se_dat),
    WM_cov = cov(d$rex - z*d$se_mod, d$rex + z*d$se_mod)[1], WM_len = cov(d$rex - z*d$se_mod, d$rex + z*d$se_mod)[2],
    WD_cov = cov(d$rex - z*d$se_dat, d$rex + z*d$se_dat)[1], WD_len = cov(d$rex - z*d$se_dat, d$rex + z*d$se_dat)[2],
    TI_cov = cov(d$ti_lo, d$ti_hi)[1], TI_len = cov(d$ti_lo, d$ti_hi)[2],
    OLD_cov = cov(d$rfirst - z*d$se_old, d$rfirst + z*d$se_old)[1]))
}
write.csv(S, file.path(OUT, "estimation_summary.csv"), row.names = FALSE); print(S, digits = 3)
