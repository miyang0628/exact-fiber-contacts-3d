## exact zeta_{1,r} and E sin(gamma) for every non-torus configuration (used by 08_analysis)
source("scripts/00_setup.R")
K <- reps(1e5, 5e3); res <- list()
for (c in c(.75, .5, .25)) res[[paste0("Esin_", c)]] <- as.list(Esin_mc(c, reps(4e6, 2e5)))
for (nm in setdiff(names(CFG), c("torus_eq","torus_mix"))) {
  cf <- CFG[[nm]]
  parts <- mclapply(seq_len(CORES), function(k) zeta_mc(cf$r, cf$lens, cf$probs, RHO, cf$box %||% c(1,1,1), cf$c %||% 1,
                    K = ceiling(K/CORES), seed = 11 + 7*k), mc.cores = CORES)
  P <- do.call(rbind, parts)
  res[[nm]] <- list(zeta = mean(P[, "zeta"]), se = sqrt(sum(P[, "se"]^2))/nrow(P), hmean = mean(P[, "hmean"]))
  cat(nm, unlist(res[[nm]]), "\n")
}
write_json(res, file.path(OUT, "zeta_configs.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
