# Paper I — R code and results

**Exact moments of fiber contact counts in three-dimensional random networks:
Steiner excluded volumes, finite-radius projections, and exact computation**
(Yang & Chun)

R (>= 4.2) with `Rcpp` and `jsonlite`; a C++ compiler is needed once to compile `R/engine.cpp`.
Run everything from this folder.

```
Rscript scripts/run_all.R                  # full study (about 15 min on 2 cores)
QUICK=1 CORES=2 Rscript scripts/run_all.R  # smoke test with small replicate counts
Rscript scripts/09_figures.R && Rscript scripts/10_paper_tables.R   # tables/figures from results/
```
`RESULTS=<dir>` redirects all inputs and outputs to another results folder.

## Contents
| Path | Content |
|---|---|
| `R/engine.cpp`, `R/engine.R` | contact kernel (closest points, gate), sweep-and-prune, cell list (Algorithm 1), brute force, torus counting, nested MC estimator of h_{1,r}; closed forms (exact q_r, lambda1/2, tip moments, zeta1 second order) |
| `R/inference.R` | exact-mean radius estimator (eq. rhat), delta-method variance |
| `scripts/01_moments.R` | Table 2 |
| `scripts/02_q_sweep.R` | Table 3, Fig. 1 |
| `scripts/03_zeta_exact.R` | Table 4, Fig. 2 |
| `scripts/04_configs.R`, `05_zeta_configs.R`, `06_analysis_configs.R` | Tables 5-6, Fig. 3 |
| `scripts/07_estimation.R` | Table 7 |
| `scripts/08_benchmark.R` | Tables 8-9, Fig. 4 (timings depend on hardware) |
| `scripts/09_figures.R` | Figures 1-4 |
| `scripts/10_paper_tables.R` | `tables/Table2...Table9_*.csv` |
| `results/` | raw replicates and summaries used in the paper |
| `tables/` | the paper's tables as CSV (paper numbering; Table 1 is a status table without data) |
| `figures/` | the paper's figures |

## Results files
- `counts_<config>.csv`: contact counts per replicate (R = 4000; torus 2000) for n = 100, ..., 3200
  (nested: N_n is the count among the first n fibers). Configurations: base (mixture {0.2,0.4},
  p1 = 0.5, r = 0.01), r005/r002/r001, equal (l = 0.3), p02/p08, three ({0.1,0.3,0.5}),
  c075/c050/c025 (oblate anisotropy), box (2 x 1 x 0.5), torus_eq, torus_mix. Margin rho = 0.01.
- `moments.csv`, `q_sweep.csv`, `zeta_exact.csv`, `zeta_configs.json`: Tables 2-4 and the exact
  zeta_{1,r} of each configuration.
- `config_analysis.csv`: per configuration and n, simulated vs exact mean and variance (and further
  distributional summaries used in Paper II).
- `estimation_replicates.csv`, `estimation_summary.csv`: radius estimators (exact, first-order,
  thin-fiber), standard errors, and interval limits (intervals are evaluated in Paper II).
- `timing.csv`, `grid.csv`, `benchmark_extra.json`: timings (R 4.3.3, Rcpp 1.0.12, Intel Xeon 2.1 GHz),
  exactness verification (240 networks, 0 mismatches), bootstrap cost.

Seeds are fixed (L'Ecuyer-CMRG streams); with the same R version and `CORES=2` a rerun reproduces the
results except timings.
