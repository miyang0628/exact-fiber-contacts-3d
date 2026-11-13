## Tables 14-15 / Fig 8: exactness check and timings of brute force, sweep-and-prune and cell list
source("scripts/00_setup.R")
r <- .01; L2 <- c(.2,.4); P2 <- c(.5,.5); set.seed(2026)
tm <- function(expr, k) {          # ms per call; each timing repeats the call until >= 50 ms elapse
  e <- substitute(expr); env <- parent.frame()
  one <- function() { reps <- 0; t0 <- proc.time()[["elapsed"]]
    repeat { eval(e, env); reps <- reps + 1; dt <- proc.time()[["elapsed"]] - t0; if (dt >= .05) break }
    1000*dt/reps }
  t <- replicate(k, one())
  c(median = median(t), q25 = quantile(t, .25, names = FALSE), q75 = quantile(t, .75, names = FALSE)) }
## verification
mism <- 0; checks <- 0
for (n in c(400L, 1600L, 3200L)) for (i in 1:reps(20, 3)) {
  f <- sample_fibers(n, L2, P2, RHO); b <- count_brute(f$X, f$T, f$L, r, c(1,1,1))
  mism <- mism + (count_sap(f$X, f$T, f$L, r, c(1,1,1))[1] != b)
  for (h in c(.05, .1, .2, .4)) { checks <- checks + 1; mism <- mism + (count_cell(f$X, f$T, f$L, r, c(1,1,1), h, -1)[1] != b) }
}
cat("verification: checks", checks, "mismatches", mism, "\n")
k <- reps(15, 3); fw <- NULL
for (n in c(400L, 800L, 1600L, 3200L, 6400L)) {
  f <- sample_fibers(n, L2, P2, RHO); s <- count_sap(f$X, f$T, f$L, r, c(1,1,1))
  fw <- rbind(fw, data.frame(regime = "fixed_window", n = n, W = 1,
    brute = if (n <= 3200) t(tm(count_brute(f$X, f$T, f$L, r, c(1,1,1)), max(3, k %/% 2))) else t(c(median=NA,q25=NA,q75=NA)),
    sap = t(tm(count_sap(f$X, f$T, f$L, r, c(1,1,1)), k)), cell = t(tm(count_cell(f$X, f$T, f$L, r, c(1,1,1), .1, -1), k)),
    candidates = s[2], contacts = s[1], pairs = n*(n-1)/2))
}
dens <- 400/Svol(L2, RHO)
for (n in c(400L, 1600L, 6400L, 25600L)) {
  W <- (n/dens)^(1/3) - max(L2) - 2*RHO; bx <- c(W, W, W)
  f <- sample_fibers(n, L2, P2, RHO, bx); s <- count_sap(f$X, f$T, f$L, r, bx)
  fw <- rbind(fw, data.frame(regime = "fixed_density", n = n, W = W,
    brute = if (n <= 6400) t(tm(count_brute(f$X, f$T, f$L, r, bx), 3)) else t(c(median=NA,q25=NA,q75=NA)),
    sap = t(tm(count_sap(f$X, f$T, f$L, r, bx), k)), cell = t(tm(count_cell(f$X, f$T, f$L, r, bx, .2, -1), k)),
    candidates = s[2], contacts = s[1], pairs = n*(n-1)/2))
}
write.csv(fw, file.path(OUT, "timing.csv"), row.names = FALSE); print(fw)
f <- sample_fibers(3200L, L2, P2, RHO); gs <- NULL
for (h in c(.025, .05, .1, .2, .4, .8)) {
  cc <- count_cell(f$X, f$T, f$L, r, c(1,1,1), h, -1); ng <- floor(2.22/h) + 1
  gs <- rbind(gs, data.frame(h = h, time = t(tm(count_cell(f$X, f$T, f$L, r, c(1,1,1), h, -1), k)),
    candidates = cc[2], registrations = cc[3], approx_mem_kB = (8*cc[3] + 16*(ng^3 + 1))/1024))
}
write.csv(gs, file.path(OUT, "grid.csv"), row.names = FALSE); print(gs)
## bootstrap cost: 200 networks at n = 1600 (generation + counting)
bc <- sapply(list(brute = function(f) count_brute(f$X, f$T, f$L, r, c(1,1,1)),
                  sap = function(f) count_sap(f$X, f$T, f$L, r, c(1,1,1)),
                  cell = function(f) count_cell(f$X, f$T, f$L, r, c(1,1,1), .1, -1)),
             function(g) { set.seed(5); system.time(for (b in 1:reps(200, 10)) g(sample_fibers(1600L, L2, P2, RHO)))[["elapsed"]] })
write_json(list(verification = list(checks = checks, mismatches = mism), bootstrap_200_networks_n1600_seconds = as.list(bc),
                R = R.version.string, Rcpp = as.character(packageVersion("Rcpp"))),
           file.path(OUT, "benchmark_extra.json"), auto_unbox = TRUE, pretty = TRUE)
print(bc)
