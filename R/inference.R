## inference.R -- exact-mean radius estimator and intervals (Section 5); base mixture design
LENS <- c(.2, .4); PROBS <- c(.5, .5); LBAR <- 0.3
SV <- Svol(LENS, RHO)
A1 <- pi*LBAR^2/SV^2; A2 <- 8*pi*LBAR/SV^2; A3 <- 32/3*pi/SV^2
qex  <- function(r) A1*r + A2*r^2 + A3*r^3                       # eq. (28)
dqex <- function(r) A1 + 2*A2*r + 3*A3*r^2
inv_ex <- function(p) { r <- p/A1; for (i in 1:50) r <- r - (qex(r) - p)/dqex(r); r }  # eq. (29)
inv_lead <- function(p) p/A1
inv_first <- function(p) { a <- A1*8/LBAR; (-A1 + sqrt(A1^2 + 4*a*p))/(2*a) }
.rg <- seq(1e-5, 0.03, length.out = 400)
.zg <- sapply(.rg, zeta1_second, lens = LENS, probs = PROBS, rho = RHO)
z2 <- function(r) approx(.rg, .zg, r, rule = 2)$y
varN <- function(r, n) { B <- n*(n-1)/2; q <- qex(r); B*(2*(n-2)*z2(r) + q*(1-q)) }
## test-inversion interval {r : |N - m(r)| <= z sd(r)}
ti_interval <- function(N, n, z) {
  B <- n*(n-1)/2
  lo <- if (N > 0) uniroot(function(r) B*qex(r) + z*sqrt(varN(r, n)) - N, c(1e-7, 0.05))$root else 0
  hi <- uniroot(function(r) B*qex(r) - z*sqrt(varN(r, n)) - N, c(1e-7, 0.05))$root
  c(lo, hi)
}
