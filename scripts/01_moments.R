## Table 3: closed-form vs Monte Carlo truncated-length and tip moments; E sin(gamma)
source("scripts/00_setup.R")
N <- reps(4e6, 2e5); set.seed(99)
res <- list()
for (cf in list(list("cube_mix", c(.2,.4), c(.5,.5), c(1,1,1)), list("cube_eq", .3, 1, c(1,1,1)),
                list("box_mix", c(.2,.4), c(.5,.5), c(2,1,.5)))) {
  nm <- cf[[1]]; lens <- cf[[2]]; probs <- cf[[3]]; box <- cf[[4]]
  f <- sample_fibers(N, lens, probs, RHO, box); l <- lC(f$X, f$T, f$L, box)
  inC <- function(P) rowSums(sweep(P, 2, 0, ">=") & sweep(P, 2, box, "<=")) == 3
  Tn <- inC(f$X + 0.5 * f$L * f$T) + inC(f$X - 0.5 * f$L * f$T)
  cf_val <- c(lam_moments(lens, probs, RHO, box), VarlC = NA, tip_moments(lens, probs, RHO, box))
  cf_val["VarlC"] <- cf_val["lambda2"] - cf_val["lambda1"]^2
  if (any(box != 1)) cf_val[4:6] <- NA
  mc <- c(mean(l), mean(l^2), var(l), mean(Tn), mean(Tn^2), mean(l * Tn))
  res[[nm]] <- data.frame(config = nm, quantity = c("lambda1","lambda2","VarlC","ET","ET2","ElT"),
                          closed_form = unname(cf_val), monte_carlo = mc)
}
es <- Esin_mc(1, N)
tab <- rbind(do.call(rbind, res), data.frame(config = "iso", quantity = "Esin", closed_form = pi/4, monte_carlo = es[1]))
write.csv(tab, file.path(OUT, "moments.csv"), row.names = FALSE); print(tab)
