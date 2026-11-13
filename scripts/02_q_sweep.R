## Table 4 / Fig 1: exact contact probability across radii, margin rho = 0.01 and rho = 0
source("scripts/00_setup.R")
lens <- c(.2,.4); probs <- c(.5,.5); n <- 3200L; R <- reps(2000, 100); B <- n*(n-1)/2
out <- NULL
for (rho in c(0.01, 0)) for (r in c(.01, .005, .002, .001)) {
  N <- unlist(par_reps(R, function() { f <- sample_fibers(n, lens, probs, rho)
        count_box(f$X, f$T, f$L, r, c(1,1,1), n) }, seed = round(1e6*r + 1e5*rho)))
  qe <- q_exact(r, lens, probs, rho); ql <- q_lead(r, lens, probs, rho)
  out <- rbind(out, data.frame(rho = rho, r = r, R = R, q_obs = mean(N)/B, se = sd(N)/sqrt(R)/B,
    q_exact = qe, q_lead = ql, obs_over_exact = mean(N)/B/qe, obs_over_lead = mean(N)/B/ql,
    first_order = 1 + 8*r/sum(lens*probs), exact_over_lead = qe/ql))
  print(tail(out, 1))
}
write.csv(out, file.path(OUT, "q_sweep.csv"), row.names = FALSE)
