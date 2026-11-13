## Tables 2-9 of Paper I as CSV files (paper numbering), from the results folder
source("scripts/00_setup.R")
TD <- file.path(ROOT, "tables"); dir.create(TD, showWarnings = FALSE)
rd <- function(f) read.csv(file.path(OUT, f)); wr <- function(x, f) write.csv(x, file.path(TD, f), row.names = FALSE)
wr(rd("moments.csv"), "Table2_truncated_moments.csv")
wr(rd("q_sweep.csv"), "Table3_contact_probability_vs_radius.csv")
z <- rd("zeta_exact.csv"); z$err_lead_pct <- 100*(z$lead/z$mc_exact - 1); z$err_mult_pct <- 100*(z$mult/z$mc_exact - 1)
z$err_second_pct <- 100*(z$second/z$mc_exact - 1); wr(z, "Table4_projection_variance.csv")
CA <- rd("config_analysis.csv")
wr(CA[CA$config != "torus_eq", c("config","n","R","nstar","mean","EN","mean_ratio","mean_ratio_se","var","VN","var_ratio","proj_frac")],
   "Table5_exact_mean_variance_design.csv")
d <- CA[CA$config %in% c("equal","torus_eq","torus_mix"), ]; B <- d$n*(d$n - 1)/2; q <- d$EN/B
wr(data.frame(config = d$config, n = d$n, sim_ratio = d$var/B/(q*(1 - q)), sim_se = d$var/B/(q*(1 - q))*sqrt(2/(d$R - 1)),
              exact_ratio = d$VN/B/(q*(1 - q))), "Table6_variance_order_cube_torus.csv")
e <- rd("estimation_summary.csv")
wr(e[, c("r","n","R","bias","bias_se","rmse","sd","first_bias","first_bias_se","first_pseudo","lead_bias","lead_pseudo","se_mod","se_dat")],
   "Table7_radius_estimators.csv")
tm <- rd("timing.csv"); wr(tm, "Table8_timing.csv"); wr(rd("grid.csv"), "Table9_cell_grid_size.csv")
cat("tables written to", TD, "\n")
