## Figures 1-4 of Paper I (exact moments) (grayscale PDFs) from the CSV results
source("scripts/00_setup.R")
FIG <- file.path(ROOT, "figures"); dir.create(FIG, showWarnings = FALSE)
rd <- function(f) read.csv(file.path(OUT, f))
P <- function(name, w, h, mf = c(1, 1)) { pdf(file.path(FIG, name), w, h, family = "serif")
  par(mfrow = mf, mar = c(3.6, 3.9, 1.2, 0.6), mgp = c(2.3, 0.7, 0), las = 1, bty = "l", cex.axis = .85, cex.lab = .95) }
ebar <- function(x, y, e, ...) arrows(x, y - e, x, y + e, angle = 90, code = 3, length = .025, ...)
## Fig 1: q ratio vs r
Q <- rd("q_sweep.csv"); P("fig_rsweep.pdf", 5, 3.4)
rr <- seq(0, .0105, length.out = 200)
plot(rr, 1 + 8*rr/.3 + 32*rr^2/(3*.09), type = "l", lwd = 1.3, xlab = "fiber radius r", ylab = expression(q[obs]/q^(0)), ylim = c(.98, 1.31))
lines(rr, 1 + 8*rr/.3, lty = 2); abline(h = 1, lty = 3, col = "grey55")
for (rho in c(.01, 0)) { d <- Q[Q$rho == rho, ]; x <- d$r + ifelse(rho == 0, 1.5e-4, 0)
  ebar(x, d$obs_over_lead, 2*d$se/d$q_lead); points(x, d$obs_over_lead, pch = ifelse(rho == 0, 22, 21), bg = ifelse(rho == 0, "white", "black"), cex = .9) }
legend("topleft", c(expression("exact " * q[r]/q^(0) * " (Thm. 4)"), expression("first order " * 1 + 8*r/bar(L)), expression("thin fiber " * q^(0)),
  expression("MC, " * rho == 0.01), expression("MC, " * rho == 0)), lty = c(1, 2, 3, NA, NA), pch = c(NA, NA, NA, 21, 22),
  pt.bg = c(NA, NA, NA, "black", "white"), col = c(1, 1, "grey55", 1, 1), bty = "n", cex = .75)
dev.off()
## Fig 2: zeta approximations
Zt <- rd("zeta_exact.csv"); P("fig_zeta.pdf", 7, 3, c(1, 2))
for (m in c("equal", "mix")) { d <- Zt[Zt$model == m, ]; band <- 200*d$mc_se/d$mc_exact
  plot(d$r, 100*(d$lead/d$mc_exact - 1), log = "x", type = "n", ylim = c(-40, 8), xlab = "r",
       ylab = expression("relative error vs. exact " * zeta[list(1, r)] * " (%)"), main = if (m == "equal") "equal length 0.3" else "mixture {0.2, 0.4}", font.main = 1, cex.main = .9)
  polygon(c(d$r, rev(d$r)), c(band, -rev(band)), col = "grey85", border = NA); abline(h = 0, col = "grey50")
  lines(d$r, 100*(d$lead/d$mc_exact - 1), type = "b", lty = 3, pch = 1)
  lines(d$r, 100*(d$mult/d$mc_exact - 1), type = "b", lty = 2, pch = 2)
  lines(d$r, 100*(d$second/d$mc_exact - 1), type = "b", lty = 1, pch = 16) }
legend("bottomleft", c("leading order", expression("multiplicative " * (1 + 8*r/bar(L))^2), "second order (Prop. 5)", "+/- 2 MC s.e."),
  lty = c(3, 2, 1, NA), pch = c(1, 2, 16, 15), col = c(1, 1, 1, "grey85"), pt.cex = c(1, 1, 1, 1.8), bty = "n", cex = .7)
dev.off()
CA <- rd("config_analysis.csv")
## Fig 3: dichotomy
P("fig_dichotomy.pdf", 5, 3.4); k <- 0
for (nm in c("equal", "torus_eq", "torus_mix")) { k <- k + 1; d <- CA[CA$config == nm, ]; B <- d$n*(d$n - 1)/2; qv <- d$EN/B
  y <- d$var/B/(qv*(1 - qv)); se <- y*sqrt(2/(d$R - 1))
  if (k == 1) plot(d$n, y, log = "xy", type = "n", ylim = c(.85, 6), xlab = "n", ylab = "variance ratio (Var / Bernoulli value)")
  lines(d$n, y, lwd = .8); ebar(d$n, y, 2*se); points(d$n, y, pch = c(21, 22, 24)[k], bg = "white") }
legend("topleft", c("cube, equal length", "torus, equal length", "torus, mixture"), pch = c(21, 22, 24), bty = "n", cex = .8); dev.off()
## Fig 4: timing, fixed density
if (file.exists(file.path(OUT, "timing.csv"))) { Tm <- rd("timing.csv"); d <- Tm[Tm$regime == "fixed_density", ]
  P("fig_celllist.pdf", 5, 3.4)
  plot(d$n, d$cell.median, log = "xy", type = "b", pch = 21, bg = "white", ylim = c(.05, 5e3), xlim = c(300, 3e4), xlab = "n", ylab = "time per network (ms)")
  lines(d$n, d$sap.median, type = "b", pch = 24, bg = "white"); ok <- !is.na(d$brute.median); lines(d$n[ok], d$brute.median[ok], type = "b", pch = 22, bg = "white")
  x <- c(300, 3e4); lines(x, d$cell.median[1]*x/400, lty = 3); lines(x, d$brute.median[1]*(x/400)^2, lty = 2)
  legend("topleft", c("brute force", "sweep-and-prune", "cell list", "slope 1 (linear)", "slope 2 (quadratic)"), pch = c(22, 24, 21, NA, NA), lty = c(1, 1, 1, 3, 2), bty = "n", cex = .75)
  dev.off() }
cat("figures written to", FIG, "\n")
cat("figures written to", FIG, "\n")
