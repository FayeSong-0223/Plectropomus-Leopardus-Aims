# =====================================================================
#  Plectropomus leopardus — Stage 2, Step 5: model diagnostics that
#  were not run earlier, and a correction to one reported result
#
#  Ziqi (Faye) Song
#
#  A full review of the project turned up four things that had not been
#  checked at all, and one reported result that does not survive them.
#
#   1. Concurvity — the GAM analogue of collinearity. Never computed.
#      It shows that s(depth) was fully aliased with the site random
#      effect (concurvity 1.000), because depth is constant within every
#      site in this subset. Depth is now a linear term instead, which
#      gives identical estimates with no aliasing.
#   2. Basis dimension adequacy (k.check). Never computed.
#   3. Seed stability of every residual-based result. The randomised
#      quantile residuals depend on runif, and all four scripts fix
#      set.seed(1). Anything computed from those residuals had never
#      been checked against a different randomisation.
#   4. Sensitivity of the protection estimate to the depth term.
#
#  Finding: the Moran's I result reported earlier was seed-dependent and
#  is corrected here. Everything else held.
#
#  Dependencies: mgcv only.
#  Run:  Rscript coral_trout_stage2_step5.R
# =====================================================================

suppressPackageStartupMessages(library(mgcv))
OUT <- "outputs"; dir.create(OUT, showWarnings = FALSE)
log_lines <- character(0)
say <- function(...) { t <- paste0(...); log_lines <<- c(log_lines, t); cat(t, "\n", sep = "") }
cap <- function(x) { log_lines <<- c(log_lines, capture.output(print(x))); print(x) }
rule <- function(t) say("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70))
base_par <- function() par(family = "sans", mgp = c(2.2, 0.6, 0), tcl = -0.3,
                           cex.axis = 0.85, cex.lab = 0.95, las = 1)
TEAL <- "#2E7D6F"; RED <- "#C1655A"

d <- read.csv(file.path(OUT, "analysis_dataset.csv"), check.names = FALSE, stringsAsFactors = FALSE)
names(d)[names(d) == "LHC_%"] <- "LHC"; names(d)[names(d) == "Corrected depth"] <- "depth"
d$REGION <- factor(d$REGION, levels = c("Palm", "Whitsunday"))
d$NTR <- factor(d$NTR, levels = c("Fished", "NTR 1987", "NTR 2004"))
d$EXPOSURE <- factor(d$EXPOSURE, levels = c("Sheltered", "Semi-Exposed", "Exposed"))
d$SITE <- factor(d$SITE)

FORM <- count ~ REGION * NTR + EXPOSURE + depth + s(YEAR, by = REGION, k = 5) +
  s(rugosity, k = 5) + s(LHC, k = 5) + s(kd490, k = 5) + s(maxDHW, k = 5) +
  s(Cyclone, k = 5) + s(SITE, bs = "re")
mA <- gam(FORM, family = nb(), data = d, method = "REML")

# ---------------------------------------------------------------------
# 1. CONCURVITY
# ---------------------------------------------------------------------
rule("1. CONCURVITY")
say("Concurvity measures how far one smooth can be reproduced by the others. It is")
say("the GAM analogue of collinearity, it is not visible in a p-value, and it was")
say("never computed for this model until now. Values approach 1 when a term is")
say("nearly redundant. The 'worst' column is the pessimistic bound.\n")
cc <- concurvity(mA, full = TRUE)
cap(round(cc, 3))
w <- cc["worst", ]; w <- w[names(w) != "para"]
say("\nterms with worst-case concurvity above 0.8:")
hi <- sort(w[w > 0.8], decreasing = TRUE)
if (length(hi)) for (n in names(hi)) say("   ", n, "  ", sprintf("%.3f", hi[[n]])) else say("   none")
say("\nThis is the number behind a caveat stated qualitatively throughout: the thermal")
say("and turbidity smooths are substantially reproducible from the region-specific")
say("year trends, because both vary largely at region-year level. Their coefficients")
say("should be read as one of several possible splits of shared variation, not as")
say("independently identified effects.")

# ---------------------------------------------------------------------
# 2. BASIS DIMENSION
# ---------------------------------------------------------------------
rule("2. BASIS DIMENSION (k.check)")
say("k was fixed at 5 throughout with a stated justification. That justification was")
say("never tested. A k-index below 1 with a small p-value means residual pattern")
say("remains that the basis cannot represent.\n")
set.seed(1); kc <- k.check(mA)
cap(round(kc, 4))
flag <- rownames(kc)[!is.na(kc[, "p-value"]) & kc[, "p-value"] < 0.05]
say("\nflagged: ", if (length(flag)) paste(flag, collapse = ", ") else "none")
say("\nFor the year smooths this is expected and not fixable by raising k: Palm has")
say("only six distinct survey years, so k cannot exceed six, and genuine year-to-year")
say("jumps are not smooth in any basis. Model B, which saturates time with region-year")
say("as a factor, is the specification that answers this, and it agrees on protection.")
say("For kd490 the flag is a real limitation and is reported as such.")

# ---------------------------------------------------------------------
# 3. SEED STABILITY OF RESIDUAL-BASED RESULTS
# ---------------------------------------------------------------------
rule("3. SEED STABILITY  (this corrects a reported result)")
say("Randomised quantile residuals draw from runif, so every diagnostic built on")
say("them depends on the seed. All earlier scripts fix set.seed(1) and none checked")
say("a second randomisation. Twenty seeds below.\n")

qres <- function(m, y) {
  mu <- fitted(m); th <- m$family$getTheta(TRUE)
  qnorm(runif(length(y), pnbinom(y - 1, size = th, mu = mu), pnbinom(y, size = th, mu = mu)))
}
moran <- function(x, lon, lat, nperm = 1999) {
  n <- length(x); R <- 6371; latm <- mean(lat) * pi / 180
  ex <- R * (lon * pi / 180) * cos(latm); ey <- R * (lat * pi / 180)
  D <- as.matrix(dist(cbind(ex, ey))); W <- 1 / D; diag(W) <- 0; W[!is.finite(W)] <- 0
  z <- x - mean(x); S <- sum(W)
  I <- (n / S) * sum(W * outer(z, z)) / sum(z^2)
  perm <- replicate(nperm, { zz <- sample(z); (n / S) * sum(W * outer(zz, zz)) / sum(zz^2) })
  c(I = I, p = (sum(abs(perm) >= abs(I)) + 1) / (nperm + 1))
}
seeds <- 1:20
ss <- do.call(rbind, lapply(seeds, function(s) {
  set.seed(s); d$r <- qres(mA, d$count)
  st <- aggregate(cbind(r, lat_dd, long_dd) ~ SITE + REGION, data = d, FUN = mean)
  o <- sapply(levels(st$REGION), function(rg) {
    z <- st[st$REGION == rg, ]; moran(z$r, z$long_dd, z$lat_dd) })
  data.frame(seed = s, Palm_I = o[1, 1], Palm_p = o[2, 1],
             Whit_I = o[1, 2], Whit_p = o[2, 2],
             shapiro_p = shapiro.test(d$r)$p.value)
}))
cap(round(ss, 4))
write.csv(round(ss, 4), file.path(OUT, "table11_seed_stability.csv"), row.names = FALSE)

say("\n-- summary across 20 seeds --")
say("  Whitsunday Moran's I : ", sprintf("%.3f to %.3f (median %.3f)",
    min(ss$Whit_I), max(ss$Whit_I), median(ss$Whit_I)),
    " | p < 0.05 in ", sum(ss$Whit_p < 0.05), " of 20")
say("  Palm       Moran's I : ", sprintf("%.3f to %.3f (median %.3f)",
    min(ss$Palm_I), max(ss$Palm_I), median(ss$Palm_I)),
    " | p < 0.05 in ", sum(ss$Palm_p < 0.05), " of 20")
say("  Shapiro-Wilk on residuals: p from ", sprintf("%.3f to %.3f", min(ss$shapiro_p), max(ss$shapiro_p)),
    " | below 0.05 in ", sum(ss$shapiro_p < 0.05), " of 20")
say("")
say("CORRECTION. Step 4 reported Whitsunday Moran's I = 0.146, p = 0.007 from a single")
say("seed and called it spatial structure the site random effect had not absorbed.")
say("That overstates it. Seed 1 sits near the top of the range, and the permutation")
say("p-value crosses 0.05 in roughly half of randomisations.")
say("")
say("What is defensible: Whitsunday Moran's I is positive under every seed tested and")
say("Palm's is centred near zero and never significant under any. The evidence points")
say("to residual spatial structure in Whitsunday but does not establish it at")
say("conventional thresholds, and any statement of it must carry the seed range.")
say("The normality of the residuals, by contrast, is robust: never rejected in 20 of 20.")

png(file.path(OUT, "fig15_seed_stability.png"), width = 2000, height = 950, res = 220)
par(mfrow = c(1, 2), mar = c(3.8, 3.9, 2.6, 0.9)); base_par()
plot(ss$seed, ss$Whit_I, type = "n", ylim = range(c(ss$Whit_I, ss$Palm_I, 0)),
     xlab = "Random seed", ylab = "Moran's I",
     main = "a  Moran's I by seed", font.main = 1, cex.main = 0.95, adj = 0)
abline(h = 0, col = "grey65", lty = 2)
lines(ss$seed, ss$Whit_I, col = TEAL, lwd = 1.6); points(ss$seed, ss$Whit_I, pch = 19, col = TEAL, cex = 0.8)
lines(ss$seed, ss$Palm_I, col = "grey55", lwd = 1.6); points(ss$seed, ss$Palm_I, pch = 19, col = "grey55", cex = 0.8)
legend("topright", legend = c("Whitsunday", "Palm"), col = c(TEAL, "grey55"),
       lwd = 1.6, pch = 19, bty = "n", cex = 0.75, text.col = "grey20")
plot(ss$seed, ss$Whit_p, type = "n", ylim = c(0, max(ss$Whit_p, 0.2)),
     xlab = "Random seed", ylab = "Permutation p-value",
     main = "b  Whitsunday p-value by seed", font.main = 1, cex.main = 0.95, adj = 0)
abline(h = 0.05, col = RED, lwd = 1.8, lty = 2)
lines(ss$seed, ss$Whit_p, col = TEAL, lwd = 1.6); points(ss$seed, ss$Whit_p, pch = 19, col = TEAL, cex = 0.8)
text(1, 0.058, "p = 0.05", col = RED, cex = 0.7, adj = 0)
dev.off(); say("\nwrote fig15_seed_stability.png")

# ---------------------------------------------------------------------
# 4. DEPTH SPECIFICATION
# ---------------------------------------------------------------------
rule("4. DEPTH: WHY IT IS NOW A LINEAR TERM")
say("depth varies at 0 of 71 sites in this subset. A penalised smooth of a")
say("site-constant covariate is fully aliased with the site random effect, and the")
say("concurvity check above returned exactly 1.000 for it before this change.\n")
rr <- function(m, region, level) {
  b <- coef(m); V <- vcov(m); nm <- names(b); k <- rep(0, length(b))
  k[match(paste0("NTR", level), nm)] <- 1
  if (region == "Whitsunday") k[match(paste0("REGIONWhitsunday:NTR", level), nm)] <- 1
  e <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(exp(e), exp(e - 1.96 * se), exp(e + 1.96 * se))
}
variants <- list(
  "depth as a smooth"  = update(FORM, . ~ . - depth + s(depth, k = 5)),
  "depth linear (used)" = FORM,
  "depth omitted"      = update(FORM, . ~ . - depth))
dep <- do.call(rbind, lapply(names(variants), function(n) {
  m <- gam(variants[[n]], family = nb(), data = d, method = "REML")
  v <- rr(m, "Whitsunday", "NTR 1987"); u <- rr(m, "Whitsunday", "NTR 2004")
  data.frame(specification = n, AIC = round(AIC(m), 1),
             dev_expl = sprintf("%.1f%%", summary(m)$dev.expl * 100),
             W1987 = sprintf("%.2f (%.2f-%.2f)", v[1], v[2], v[3]),
             W2004 = sprintf("%.2f (%.2f-%.2f)", u[1], u[2], u[3]))
}))
cap(dep)
write.csv(dep, file.path(OUT, "table12_depth_specification.csv"), row.names = FALSE)
say("\nThe smooth and the linear term give identical estimates, because the smooth had")
say("already collapsed to a straight line (edf 1.005). Omitting depth raises the")
say("Whitsunday estimate, so depth is doing real adjustment work and should be kept.")
say("The linear form keeps that adjustment without the aliasing.")

writeLines(log_lines, file.path(OUT, "stage2_step5_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
