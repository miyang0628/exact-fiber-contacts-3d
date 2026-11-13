## Table 5 / Fig 2: exact zeta_{1,r} by nested MC vs leading, multiplicative, second order
source("scripts/00_setup.R")
K <- reps(2e5, 1e4); out <- NULL
for (m in list(list("equal", .3, 1), list("mix", c(.2,.4), c(.5,.5)))) for (r in c(.001,.002,.005,.01)) {
  lens <- m[[2]]; probs <- m[[3]]; Lb <- sum(lens*probs)
  ## split K over cores (independent seeds), combine unbiased pieces
  parts <- mclapply(seq_len(CORES), function(k) zeta_mc(r, lens, probs, RHO, K = ceiling(K/CORES), seed = 1000*k + round(1e5*r)),
                    mc.cores = CORES)
  P <- do.call(rbind, parts); z <- mean(P[, "zeta"]); se <- sqrt(sum(P[, "se"]^2))/nrow(P)
  zl <- zeta1_lead(r, lens, probs, RHO)
  out <- rbind(out, data.frame(model = m[[1]], r = r, lead = zl, mult = zl*(1+8*r/Lb)^2,
    second = zeta1_second(r, lens, probs, RHO), mc_exact = z, mc_se = se, q_exact = q_exact(r, lens, probs, RHO)))
  print(tail(out, 1))
}
write.csv(out, file.path(OUT, "zeta_exact.csv"), row.names = FALSE)
