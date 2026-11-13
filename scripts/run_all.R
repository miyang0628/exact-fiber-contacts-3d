## Run every experiment in order (from the project root):  Rscript scripts/run_all.R
## Quick smoke test:  QUICK=1 Rscript scripts/run_all.R      Cores:  CORES=8 Rscript ...
for (s in c("01_moments","02_q_sweep","03_zeta_exact","04_configs","05_zeta_configs","06_analysis_configs","07_estimation","08_benchmark","09_figures","10_paper_tables")) {
  cat("\n=====", s, "=====\n"); t0 <- Sys.time()
  system2("Rscript", file.path("scripts", paste0(s, ".R")))
  cat("elapsed:", format(Sys.time() - t0), "\n")
}
