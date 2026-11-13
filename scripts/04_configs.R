## Tables 6-8, 12 / Figs 3-7: nested counts N_n (n in NG) for every configuration
## Each replicate draws 3200 fibers; N_n is the count among the first n (prefix sample).
source("scripts/00_setup.R")
args <- commandArgs(TRUE); which <- if (length(args)) strsplit(args[1], ",")[[1]] else names(CFG)
for (nm in which) {
  cf <- CFG[[nm]]; box <- cf$box %||% c(1,1,1); c <- cf$c %||% 1; tor <- isTRUE(cf$torus)
  R <- if (tor) reps(2000, 40) else reps(4000, 100)
  t0 <- Sys.time()
  A <- do.call(rbind, par_reps(R, function() {
    f <- sample_fibers(max(NG), cf$lens, cf$probs, RHO, box, c, tor)
    if (tor) count_torus(f$X, f$T, f$L, cf$r, NG) else count_box(f$X, f$T, f$L, cf$r, box, NG)
  }, seed = sum(utf8ToInt(nm)) * 101))
  colnames(A) <- paste0("n", NG)
  write.csv(A, file.path(OUT, sprintf("counts_%s.csv", nm)), row.names = FALSE)
  cat(nm, nrow(A), round(colMeans(A), 2), format(Sys.time() - t0), "\n")
}
