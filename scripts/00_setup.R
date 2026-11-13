## common setup for all experiment scripts
ROOT <- normalizePath(Sys.getenv("FIBER_ROOT", "."))
if (!file.exists(file.path(ROOT, "R", "engine.R"))) ROOT <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path(ROOT, "R", "engine.R"))
suppressPackageStartupMessages({ library(parallel); library(jsonlite) })
RNGkind("L'Ecuyer-CMRG")
CORES <- as.integer(Sys.getenv("CORES", max(1, detectCores() - 0)))
QUICK <- Sys.getenv("QUICK", "0") == "1"      # QUICK=1: small replicate counts for testing
OUT <- Sys.getenv("RESULTS", file.path(ROOT, "results")); dir.create(OUT, showWarnings = FALSE)
reps <- function(full, quick) if (QUICK) quick else full
## parallel replicate helper with reproducible L'Ecuyer streams
par_reps <- function(R, FUN, seed) {
  set.seed(seed); chunks <- split(seq_len(R), rep_len(seq_len(CORES), R))
  res <- mclapply(chunks, function(ix) lapply(ix, function(i) FUN()), mc.cores = CORES, mc.set.seed = TRUE)
  unlist(res, recursive = FALSE)
}
