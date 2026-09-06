# =====================================================================
#  Plectropomus leopardus — Stage 2, Step 1: the simplest defensible model
#  Palm and Whitsunday inshore reefs, Great Barrier Reef, 2007-2018
#
#  Ziqi (Faye) Song
#
#  Fits a negative binomial mixed model with a site-level random intercept:
#
#      count ~ REGION * NTR + EXPOSURE + factor(YEAR) + s(SITE, bs = "re")
#
#  No habitat or environmental covariates yet, and no smooth terms. The
#  purpose of this step is to establish that the model machinery is sound
#  and the diagnostics are clean BEFORE anything interesting goes in, and
#  to see whether the protection x region pattern in Figure 2 survives
#  adjustment for wave exposure.
#
#  Dependencies: mgcv only, which ships with every R installation.
#  gam() with bs = "re" fits the random effect by REML, so this is a
#  mixed model in the usual sense rather than a smoothing device.
#
#  Run:  Rscript coral_trout_stage2_model.R    (from the project directory)
#  Requires: outputs/analysis_dataset.csv, produced by coral_trout_stage1.R
# =====================================================================

set.seed(1)
suppressPackageStartupMessages(library(mgcv))
DATA <- "data"; OUT <- "outputs"
dir.create(OUT, showWarnings = FALSE)

log_lines <- character(0)
say <- function(...) { t <- paste0(...); log_lines <<- c(log_lines, t); cat(t, "\n", sep = "") }
cap <- function(x) { log_lines <<- c(log_lines, capture.output(print(x))); print(x) }
rule <- function(t) say("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70))

COL_NTR <- c("Fished" = "#C1655A", "NTR 1987" = "#2E7D6F", "NTR 2004" = "#4A6FA5")
base_par <- function() par(family = "sans", mgp = c(2.2, 0.6, 0), tcl = -0.3,
                           cex.axis = 0.85, cex.lab = 0.95, las = 1)

# ---------------------------------------------------------------------
# 1. LOAD THE PREPARED DATA
# ---------------------------------------------------------------------
rule("1. LOAD")
f <- file.path(OUT, "analysis_dataset.csv")
if (!file.exists(f))
  stop("outputs/analysis_dataset.csv not found. Run coral_trout_stage1.R first.")
d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

d$REGION   <- factor(d$REGION,   levels = c("Palm", "Whitsunday"))
d$NTR      <- factor(d$NTR,      levels = c("Fished", "NTR 1987", "NTR 2004"))
d$EXPOSURE <- factor(d$EXPOSURE, levels = c("Sheltered", "Semi-Exposed", "Exposed"))
d$SITE     <- factor(d$SITE)
d$YEARF    <- factor(d$YEAR)
d$RY       <- factor(paste(d$REGION, d$YEAR))

say("rows: ", nrow(d), "   sites: ", nlevels(d$SITE), "   region-years: ", nlevels(d$RY))
say("response is the recovered integer count; survey area is constant, so no offset")
say("counts range ", min(d$count), " to ", max(d$count), ", mean ", sprintf("%.2f", mean(d$count)),
    ", variance ", sprintf("%.2f", var(d$count)))
say("variance / mean = ", sprintf("%.2f", var(d$count) / mean(d$count)),
    "  -> overdispersed relative to Poisson, so negative binomial rather than Poisson")

# ---------------------------------------------------------------------
# 2. FIT
# ---------------------------------------------------------------------
rule("2. FIT")
say("Model 1 (primary):   count ~ REGION*NTR + EXPOSURE + YEAR + s(SITE, re)")
say("Model 2 (sensitivity): the same without EXPOSURE\n")
say("The comparison is the point. Protection is confounded with exposure, and")
say("the confounding is far worse in Whitsunday than in Palm (Stage 1, section 3b).")
say("If the protection estimate moves a lot between these two fits, the unadjusted")
say("contrast in Figure 2 was substantially an exposure contrast.\n")

m1 <- gam(count ~ REGION * NTR + EXPOSURE + YEARF + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")
m2 <- gam(count ~ REGION * NTR + YEARF + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")

theta <- m1$family$getTheta(TRUE)
say("negative binomial theta (estimated): ", sprintf("%.3f", theta))
say("  small theta = strong overdispersion; theta -> Inf recovers Poisson")
say("deviance explained, model 1: ", sprintf("%.1f%%", summary(m1)$dev.expl * 100))
say("deviance explained, model 2: ", sprintf("%.1f%%", summary(m2)$dev.expl * 100))

say("\n-- random effect --")
vc <- gam.vcomp(m1, rescale = FALSE)
cap(vc)
say("\nThe site standard deviation is on the log scale. A value of s means that")
say("a site one standard deviation above average holds exp(s) times the density")
say("of an average site, after everything else in the model is accounted for.")

say("\n-- parametric coefficients, model 1 --")
cap(round(summary(m1)$p.table, 4))

# ---------------------------------------------------------------------
# 3. PROTECTION EFFECTS AS RATE RATIOS
# ---------------------------------------------------------------------
rule("3. PROTECTION EFFECTS")
say("Reported as rate ratios against Fished sites within the same region.")
say("A ratio of 2 means twice the density. The interval is 95% Wald.\n")

# build contrast vectors against the fitted coefficient names
rr <- function(m, region, level) {
  b <- coef(m); V <- vcov(m); k <- rep(0, length(b)); nm <- names(b)
  main <- paste0("NTR", level)
  if (!(main %in% nm)) stop("coefficient not found: ", main)
  k[match(main, nm)] <- 1
  if (region == "Whitsunday") {
    ix <- paste0("REGIONWhitsunday:NTR", level)
    if (ix %in% nm) k[match(ix, nm)] <- 1
  }
  est <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(ratio = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

grid <- expand.grid(region = c("Palm", "Whitsunday"),
                    level  = c("NTR 1987", "NTR 2004"), stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  a <- rr(m1, grid$region[i], grid$level[i])
  b <- rr(m2, grid$region[i], grid$level[i])
  data.frame(region = grid$region[i], protection = grid$level[i],
             adj_ratio = a[["ratio"]], adj_lo = a[["lo"]], adj_hi = a[["hi"]],
             unadj_ratio = b[["ratio"]], unadj_lo = b[["lo"]], unadj_hi = b[["hi"]])
}))
res[, 3:8] <- round(res[, 3:8], 2)
cap(res)
write.csv(res, file.path(OUT, "table3_protection_rate_ratios.csv"), row.names = FALSE)

say("\nadj   = adjusted for wave exposure (model 1)")
say("unadj = not adjusted for exposure   (model 2)")

# ---------------------------------------------------------------------
# 4. DIAGNOSTICS  (randomised quantile residuals, Dunn & Smyth 1996)
# ---------------------------------------------------------------------
rule("4. DIAGNOSTICS")
say("Deviance residuals are unreliable for discrete responses with small counts,")
say("so randomised quantile residuals are used instead: under a correct model")
say("they are exactly standard normal, which makes a QQ plot interpretable.\n")

qres <- function(m, y) {
  mu <- fitted(m); th <- m$family$getTheta(TRUE)
  a <- pnbinom(y - 1, size = th, mu = mu)
  b <- pnbinom(y,     size = th, mu = mu)
  qnorm(runif(length(y), a, b))
}
r1 <- qres(m1, d$count)
say("randomised quantile residuals: mean ", sprintf("%+.3f", mean(r1)),
    "  sd ", sprintf("%.3f", sd(r1)))
sw <- shapiro.test(r1)
say("Shapiro-Wilk on the residuals: W = ", sprintf("%.4f", sw$statistic),
    ", p = ", sprintf("%.4f", sw$p.value))
say("  (a formal test is over-powered at n = ", length(r1),
    "; the QQ plot in Figure 4 is the thing to judge)")

ry_mean <- tapply(r1, d$RY, mean)
say("\nlargest mean residual by region-year (structure the model is missing):")
cap(round(sort(ry_mean)[c(1:3, (length(ry_mean) - 2):length(ry_mean))], 3))

site_mean <- tapply(r1, d$SITE, mean)
say("\nsd of site-mean residuals: ", sprintf("%.3f", sd(site_mean)),
    "  (near 1/sqrt(obs per site) if the site random effect is doing its job)")
say("mean observations per site: ", sprintf("%.1f", nrow(d) / nlevels(d$SITE)),
    "  -> expected sd about ", sprintf("%.3f", 1 / sqrt(nrow(d) / nlevels(d$SITE))))

## ---- Figure 4: diagnostics -------------------------------------------
png(file.path(OUT, "fig4_model_diagnostics.png"), width = 2100, height = 1650, res = 230)
par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.3, 0.8)); base_par()

qqnorm(r1, main = "a  Quantile residuals", font.main = 1, cex.main = 1.0, adj = 0,
       pch = 16, col = adjustcolor("grey25", 0.45), cex = 0.6,
       xlab = "Theoretical quantiles", ylab = "Observed quantiles")
qqline(r1, col = "#C1655A", lwd = 2)

plot(log(fitted(m1)), r1, pch = 16, col = adjustcolor("grey25", 0.45), cex = 0.6,
     xlab = "Fitted value (log scale)", ylab = "Quantile residual",
     main = "b  Residuals vs fitted", font.main = 1, cex.main = 1.0, adj = 0)
abline(h = 0, col = "#C1655A", lwd = 1.6)
lines(lowess(log(fitted(m1)), r1), col = "#2E7D6F", lwd = 2)

boxplot(r1 ~ d$RY, las = 2, cex.axis = 0.5, col = "grey93", border = "grey40",
        xlab = "", ylab = "Quantile residual",
        main = "c  Residuals by region-year", font.main = 1, cex.main = 1.0, adj = 0)
abline(h = 0, col = "#C1655A", lwd = 1.6)

hist(site_mean, breaks = 18, col = "grey90", border = "white",
     xlab = "Mean residual per site", ylab = "Number of sites",
     main = "d  Site-level residual means", font.main = 1, cex.main = 1.0, adj = 0)
abline(v = 0, col = "#C1655A", lwd = 1.6)
dev.off(); say("\nwrote fig4_model_diagnostics.png")

## ---- Figure 5: protection effects ------------------------------------
png(file.path(OUT, "fig5_protection_effects.png"), width = 2000, height = 1150, res = 230)
par(mar = c(3.8, 8.6, 2.6, 1.2)); base_par()
res$lab <- paste0(res$region, ", ", res$protection)
ord <- rev(seq_len(nrow(res)))
xr <- range(c(res$adj_lo, res$adj_hi, res$unadj_lo, res$unadj_hi, 1))
plot(NA, xlim = xr, ylim = c(0.5, nrow(res) + 0.5), yaxt = "n", log = "x",
     xlab = "Density ratio vs fished sites in the same region", ylab = "",
     main = "Protection effect, with and without adjustment for wave exposure",
     font.main = 1, cex.main = 0.95, adj = 0)
axis(2, at = ord, labels = res$lab, cex.axis = 0.8)
abline(v = 1, col = "grey55", lty = 2)
grid(nx = NULL, ny = NA, col = "grey93", lty = 1)
off <- 0.16
segments(res$unadj_lo, ord + off, res$unadj_hi, ord + off, col = "grey60", lwd = 2)
points(res$unadj_ratio, ord + off, pch = 21, bg = "white", col = "grey45", cex = 1.1)
segments(res$adj_lo, ord - off, res$adj_hi, ord - off, col = "#2E7D6F", lwd = 2.4)
points(res$adj_ratio, ord - off, pch = 19, col = "#2E7D6F", cex = 1.2)
legend("topleft", legend = c("adjusted for exposure", "unadjusted"),
       col = c("#2E7D6F", "grey55"), pch = c(19, 21), pt.bg = "white",
       lwd = c(2.4, 2), bty = "n", cex = 0.78, text.col = "grey20")
dev.off(); say("wrote fig5_protection_effects.png")

# ---------------------------------------------------------------------
# 5. WHAT THIS STEP DOES AND DOES NOT SETTLE
# ---------------------------------------------------------------------
rule("5. READING THIS STEP")
say("Settled: whether the model machinery is adequate — distribution, random")
say("effect, residual behaviour — and how far the protection contrast moves when")
say("wave exposure is adjusted for.")
say("")
say("Not settled: anything about environmental drivers, which are not in the model")
say("yet; and anything causal. Protection was not assigned at random, and adjusting")
say("for exposure removes one known confounder rather than all of them.")
say("")
say("Decision before Step 2: habitat covariates are added next. If reserves were")
say("PLACED in different habitat, adjusting for habitat is correct. If protection")
say("CHANGED the habitat, habitat is a mediator and adjusting for it would remove")
say("part of the effect being measured. For coral trout over this period the first")
say("is far more plausible, so Step 2 adjusts. Recorded here so the choice is not")
say("revisited after seeing which way the coefficient moves.")

writeLines(log_lines, file.path(OUT, "stage2_step1_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
