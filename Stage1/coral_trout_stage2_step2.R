# =====================================================================
#  Plectropomus leopardus — Stage 2, Step 2: habitat, environment,
#  and where the variation actually sits
#  Palm and Whitsunday inshore reefs, Great Barrier Reef, 2007-2018
#
#  Ziqi (Faye) Song
#
#  Step 1 established that a negative binomial mixed model with a site
#  random intercept is adequate, and that the regional difference in the
#  protection effect survives adjustment for wave exposure. Two things
#  follow from its diagnostics:
#
#   (a) residuals by region-year showed clear structure, so the year
#       effect must be allowed to differ between regions; and
#   (b) habitat and environmental covariates now go in.
#
#  These two requirements conflict. Region-year fixed effects would give
#  the best possible adjustment for time, but maxDHW and SSTmean carry
#  99% and 98% of their variance between region-years (Stage 1, table 2),
#  so such a model cannot estimate a thermal effect at all. Rather than
#  choose silently, two models are fitted:
#
#   Model A (drivers)    region-specific smooth year trends, leaving
#                        year-to-year deviation available to the
#                        environmental covariates.
#   Model B (time-saturated)  region-year as a fixed factor. The cleanest
#                        possible protection estimate; thermal covariates
#                        deliberately excluded because they are not
#                        identifiable here.
#
#  Environmental effects are read from A only, protection from both, and
#  the difference between them is reported rather than hidden.
#
#  Dependencies: mgcv only (ships with R).
#  Run:  Rscript coral_trout_stage2_step2.R
#  Requires: outputs/analysis_dataset.csv from coral_trout_stage1.R
# =====================================================================

set.seed(1)
suppressPackageStartupMessages(library(mgcv))
OUT <- "outputs"; dir.create(OUT, showWarnings = FALSE)

log_lines <- character(0)
say <- function(...) { t <- paste0(...); log_lines <<- c(log_lines, t); cat(t, "\n", sep = "") }
cap <- function(x) { log_lines <<- c(log_lines, capture.output(print(x))); print(x) }
rule <- function(t) say("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70))
base_par <- function() par(family = "sans", mgp = c(2.2, 0.6, 0), tcl = -0.3,
                           cex.axis = 0.85, cex.lab = 0.95, las = 1)
TEAL <- "#2E7D6F"; RED <- "#C1655A"; BLUE <- "#4A6FA5"; AMBER <- "#B4813A"

# ---------------------------------------------------------------------
# 1. LOAD
# ---------------------------------------------------------------------
rule("1. LOAD")
f <- file.path(OUT, "analysis_dataset.csv")
if (!file.exists(f)) stop("outputs/analysis_dataset.csv not found. Run coral_trout_stage1.R first.")
d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
names(d)[names(d) == "LHC_%"]          <- "LHC"
names(d)[names(d) == "LCC_%"]          <- "LCC"
names(d)[names(d) == "Corrected depth"] <- "depth"

d$REGION   <- factor(d$REGION,   levels = c("Palm", "Whitsunday"))
d$NTR      <- factor(d$NTR,      levels = c("Fished", "NTR 1987", "NTR 2004"))
d$EXPOSURE <- factor(d$EXPOSURE, levels = c("Sheltered", "Semi-Exposed", "Exposed"))
d$SITE     <- factor(d$SITE)
d$RY       <- factor(paste(d$REGION, d$YEAR))
say("rows ", nrow(d), " | sites ", nlevels(d$SITE), " | region-years ", nlevels(d$RY))
say("covariates carried forward (chosen in Stage 1, not by selection):")
say("  habitat      rugosity, LHC (live hard coral), depth")
say("  environment  kd490 (water clarity), maxDHW (thermal), Cyclone")
say("  management   NTR x REGION, adjusted for EXPOSURE")

# ---------------------------------------------------------------------
# 2. FIT
# ---------------------------------------------------------------------
rule("2. FIT")
say("Smooth basis dimension is held to k = 5 throughout. Palm has only six survey")
say("years, so anything larger invites a wiggly fit the data cannot support.")
say("Effective degrees of freedom are reported so over-fitting is visible.\n")

mA <- gam(count ~ REGION * NTR + EXPOSURE +
            s(YEAR, by = REGION, k = 5) +
            s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
            s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) +
            s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")

mB <- gam(count ~ REGION * NTR + EXPOSURE + RY +
            s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
            s(kd490, k = 5) + s(Cyclone, k = 5) +
            s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")

for (nm in c("A", "B")) {
  m <- get(paste0("m", nm)); s <- summary(m)
  say("Model ", nm, ": theta ", sprintf("%.2f", m$family$getTheta(TRUE)),
      " | deviance explained ", sprintf("%.1f%%", s$dev.expl * 100),
      " | REML ", sprintf("%.1f", m$gcv.ubre), " | n = ", nrow(d))
}

say("\n-- Model A, smooth terms --")
cap(round(summary(mA)$s.table, 4))
say("\nedf near 1 means the fitted relationship is effectively a straight line.")
say("A smooth whose edf collapses to 1 and whose p-value is large is doing nothing.")

say("\n-- variance components --")
cap(gam.vcomp(mA, rescale = FALSE))

# ---------------------------------------------------------------------
# 3. PROTECTION ACROSS SPECIFICATIONS
# ---------------------------------------------------------------------
rule("3. IS THE PROTECTION EFFECT STABLE?")
say("The estimate is worth little if it moves with every change of specification.")
say("Four models are compared, from the Step 1 baseline to the fully adjusted fits.\n")

m0 <- gam(count ~ REGION * NTR + EXPOSURE + factor(YEAR) + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")
mH <- gam(count ~ REGION * NTR + EXPOSURE + factor(YEAR) +
            s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")

rr <- function(m, region, level) {
  b <- coef(m); V <- vcov(m); nm <- names(b); k <- rep(0, length(b))
  k[match(paste0("NTR", level), nm)] <- 1
  ix <- paste0("REGIONWhitsunday:NTR", level)
  if (region == "Whitsunday" && ix %in% nm) k[match(ix, nm)] <- 1
  e <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(exp(e), exp(e - 1.96 * se), exp(e + 1.96 * se))
}
mods <- list("1 baseline (Step 1)" = m0, "2 + habitat" = mH,
             "3 + environment, smooth year (A)" = mA, "4 region-year saturated (B)" = mB)
grid <- expand.grid(region = c("Palm", "Whitsunday"),
                    level = c("NTR 1987", "NTR 2004"), stringsAsFactors = FALSE)
stab <- do.call(rbind, lapply(names(mods), function(nm)
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    v <- rr(mods[[nm]], grid$region[i], grid$level[i])
    data.frame(model = nm, region = grid$region[i], protection = grid$level[i],
               ratio = v[1], lo = v[2], hi = v[3])
  }))))
stab[, 4:6] <- round(stab[, 4:6], 2)
cap(stab[order(stab$region, stab$protection, stab$model), ])
write.csv(stab, file.path(OUT, "table4_protection_stability.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4. VARIANCE DECOMPOSITION ON THE LINK SCALE
# ---------------------------------------------------------------------
rule("4. WHERE THE VARIATION SITS")
say("Each block contributes a vector of values to the linear predictor. The variance")
say("of each block is its unique contribution; the pairwise covariances are what the")
say("blocks share. Variances and covariances together sum exactly to the variance of")
say("the linear predictor, so nothing is allocated arbitrarily and the answer does")
say("not depend on the order terms were entered.\n")

tm <- predict(mA, type = "terms")
cn <- colnames(tm)
blk <- list(
  Management  = grep("NTR", cn, value = TRUE),
  Habitat     = grep("s\\((rugosity|LHC|depth)\\)", cn, value = TRUE),
  Environment = grep("s\\((kd490|maxDHW|Cyclone)\\)", cn, value = TRUE),
  `Space/time`= c(grep("^REGION$|s\\(YEAR\\)", cn, value = TRUE),
                  grep("YEAR", cn, value = TRUE), grep("^EXPOSURE$", cn, value = TRUE)),
  `Site (RE)` = grep("s\\(SITE\\)", cn, value = TRUE))
blk <- lapply(blk, unique)
blk <- lapply(blk, function(x) x[x %in% cn])
say("term-to-block assignment:")
for (b in names(blk)) say("  ", b, ": ", paste(blk[[b]], collapse = ", "))
used <- unlist(blk); miss <- setdiff(cn, used)
if (length(miss)) say("  UNASSIGNED (check): ", paste(miss, collapse = ", "))

B <- sapply(blk, function(k) if (length(k)) rowSums(tm[, k, drop = FALSE]) else rep(0, nrow(tm)))
S <- cov(B); tot <- sum(S)
uni <- diag(S) / tot
say("\n-- unique contribution of each block (share of linear-predictor variance) --")
cap(round(sort(uni, decreasing = TRUE), 3))
say("\n-- shared components (pairwise covariance / total) --")
sh <- S; diag(sh) <- NA
shd <- as.data.frame(as.table(round(sh / tot, 3)))
shd <- shd[!is.na(shd$Freq) & as.character(shd$Var1) < as.character(shd$Var2), ]
shd <- shd[order(-abs(shd$Freq)), ]; names(shd) <- c("block_1", "block_2", "shared_share")
cap(head(shd, 6))
say("\nsum of unique shares: ", sprintf("%.3f", sum(uni)),
    "   sum of 2x shared: ", sprintf("%.3f", 1 - sum(uni)))
say("Positive shared values mean two blocks move together and observation alone")
say("cannot separate their contributions. Negative values mean they offset.")
vd <- data.frame(block = names(uni), unique_share = round(as.numeric(uni), 4))
write.csv(vd, file.path(OUT, "table5_variance_decomposition.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 5. FIGURES
# ---------------------------------------------------------------------
rule("5. FIGURES")

## Figure 6 — partial effects
png(file.path(OUT, "fig6_partial_effects.png"), width = 2100, height = 1400, res = 220)
par(mfrow = c(2, 3), mar = c(3.6, 3.8, 2.4, 0.8)); base_par()
sm <- c("s(rugosity)", "s(LHC)", "s(depth)", "s(kd490)", "s(maxDHW)", "s(Cyclone)")
lab <- c("Rugosity index", "Live hard coral (%)", "Depth (m)",
         "kd490 (turbidity)", "Max degree heating weeks", "Cyclone exposure index")
st <- summary(mA)$s.table
for (i in seq_along(sm)) {
  j <- match(sm[i], rownames(st))
  plot(mA, select = which(sapply(mA$smooth, function(z) z$label) == sm[i]),
       shade = TRUE, shade.col = adjustcolor(TEAL, 0.20), col = TEAL, lwd = 2,
       xlab = lab[i], ylab = "Effect on log density", rug = TRUE,
       main = sprintf("%s   edf %.2f, p %s", letters[i],
                      st[j, "edf"],
                      ifelse(st[j, "p-value"] < 0.001, "< 0.001",
                             sprintf("= %.3f", st[j, "p-value"]))),
       font.main = 1, cex.main = 0.95, adj = 0)
  abline(h = 0, col = "grey60", lty = 2)
}
dev.off(); say("wrote fig6_partial_effects.png")

## Figure 7 — variance decomposition
png(file.path(OUT, "fig7_variance_decomposition.png"), width = 1900, height = 1150, res = 220)
par(mar = c(3.8, 8.0, 2.6, 1.0)); base_par()
o <- sort(uni)
bp <- barplot(o, horiz = TRUE, col = TEAL, border = NA, xlim = c(0, max(o) * 1.25),
              xlab = "Share of linear-predictor variance (unique)",
              main = "Where the explained variation sits", font.main = 1,
              cex.main = 0.98, adj = 0, las = 1, cex.names = 0.82)
text(o, bp, sprintf(" %.3f", o), pos = 4, cex = 0.8, col = "grey25", xpd = NA)
dev.off(); say("wrote fig7_variance_decomposition.png")

## Figure 8 — protection stability across specifications
png(file.path(OUT, "fig8_protection_stability.png"), width = 2000, height = 1250, res = 220)
par(mfrow = c(1, 2), mar = c(3.9, 3.9, 2.6, 0.8)); base_par()
cols <- c(RED, AMBER, TEAL, BLUE)
for (rg in c("Palm", "Whitsunday")) {
  s <- stab[stab$region == rg, ]
  s$grp <- paste(s$protection)
  xs <- seq_len(length(unique(s$model)))
  plot(NA, xlim = c(0.6, length(xs) + 0.4), ylim = range(c(stab$lo, stab$hi, 1)), log = "y",
       xaxt = "n", xlab = "", ylab = "Density ratio vs fished",
       main = rg, font.main = 1, cex.main = 1.0, adj = 0)
  axis(1, at = xs, labels = c("base", "+hab", "+env", "sat"), cex.axis = 0.8)
  abline(h = 1, col = "grey55", lty = 2); grid(nx = NA, ny = NULL, col = "grey93")
  for (k in seq_along(unique(s$protection))) {
    lv <- unique(s$protection)[k]; z <- s[s$protection == lv, ]
    z <- z[match(names(mods), z$model), ]
    off <- (k - 1.5) * 0.13
    segments(xs + off, z$lo, xs + off, z$hi, col = c(TEAL, BLUE)[k], lwd = 2)
    lines(xs + off, z$ratio, col = c(TEAL, BLUE)[k], lwd = 1.6, lty = 3)
    points(xs + off, z$ratio, pch = 19, col = c(TEAL, BLUE)[k], cex = 1.05)
  }
  if (rg == "Palm") legend("topleft", legend = unique(s$protection), col = c(TEAL, BLUE),
                           pch = 19, lwd = 2, bty = "n", cex = 0.78, text.col = "grey20")
}
dev.off(); say("wrote fig8_protection_stability.png")

# ---------------------------------------------------------------------
# 6. DIAGNOSTICS FOR MODEL A
# ---------------------------------------------------------------------
rule("6. DIAGNOSTICS")
qres <- function(m, y) {
  mu <- fitted(m); th <- m$family$getTheta(TRUE)
  qnorm(runif(length(y), pnbinom(y - 1, size = th, mu = mu), pnbinom(y, size = th, mu = mu)))
}
rA <- qres(mA, d$count)
say("Model A quantile residuals: mean ", sprintf("%+.3f", mean(rA)),
    " sd ", sprintf("%.3f", sd(rA)),
    " | Shapiro-Wilk p = ", sprintf("%.3f", shapiro.test(rA)$p.value))
ry <- tapply(rA, d$RY, mean)
say("region-year mean residuals now range ", sprintf("%+.3f", min(ry)),
    " to ", sprintf("%+.3f", max(ry)),
    "  (Step 1 range was -0.652 to +0.503)")
say("sd of region-year mean residuals: ", sprintf("%.3f", sd(ry)))

# ---------------------------------------------------------------------
# 7. IS THE THERMAL EFFECT REAL, OR IS IT WHITSUNDAY 2017?
# ---------------------------------------------------------------------
rule("7. THE THERMAL SMOOTH, INTERROGATED")
say("s(maxDHW) is the strongest environmental term, but the fitted curve RISES")
say("above about 4 degree heating weeks, which would mean heat stress increases")
say("coral trout density. That is implausible, and the rug plot in Figure 6e shows")
say("a gap in the data between roughly 4 and 5.8 DHW. The upper limb of the curve")
say("therefore rests on very few region-years.\n")

dhw_ry <- tapply(d$maxDHW, d$RY, mean)
say("region-years with maxDHW above 4:")
cap(round(dhw_ry[dhw_ry > 4], 2))
say("\nn observations above 4 DHW: ", sum(d$maxDHW > 4), " of ", nrow(d))

d2 <- droplevels(d[!(d$REGION == "Whitsunday" & d$YEAR == 2017), ])
mA2 <- gam(count ~ REGION * NTR + EXPOSURE +
             s(YEAR, by = REGION, k = 5) +
             s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
             s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) +
             s(SITE, bs = "re"),
           family = nb(), data = d2, method = "REML")
st2 <- summary(mA2)$s.table
say("\n-- Model A refitted WITHOUT the Whitsunday 2017 survey (n = ", nrow(d2), ") --")
cap(round(st2, 4))

j1 <- match("s(maxDHW)", rownames(summary(mA)$s.table))
j2 <- match("s(maxDHW)", rownames(st2))
say("\ns(maxDHW) with 2017    : edf ", sprintf("%.2f", summary(mA)$s.table[j1, "edf"]),
    ", p ", sprintf("%.4f", summary(mA)$s.table[j1, "p-value"]))
say("s(maxDHW) without 2017 : edf ", sprintf("%.2f", st2[j2, "edf"]),
    ", p ", sprintf("%.4f", st2[j2, "p-value"]))

png(file.path(OUT, "fig9_thermal_sensitivity.png"), width = 1900, height = 950, res = 220)
par(mfrow = c(1, 2), mar = c(3.7, 3.9, 2.5, 0.8)); base_par()
for (mm in list(list(m = mA,  t = "a  All data"),
                list(m = mA2, t = "b  Whitsunday 2017 removed"))) {
  k <- which(sapply(mm$m$smooth, function(z) z$label) == "s(maxDHW)")
  plot(mm$m, select = k, shade = TRUE, shade.col = adjustcolor(TEAL, 0.20),
       col = TEAL, lwd = 2, rug = TRUE, ylim = c(-1.2, 1.2),
       xlab = "Max degree heating weeks", ylab = "Effect on log density",
       main = mm$t, font.main = 1, cex.main = 1.0, adj = 0)
  abline(h = 0, col = "grey60", lty = 2)
}
dev.off(); say("\nwrote fig9_thermal_sensitivity.png")

say("\nRead this comparison before quoting any thermal effect. If the curve changes")
say("shape when one survey is removed, the effect is that survey, not the covariate.")

writeLines(log_lines, file.path(OUT, "stage2_step2_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
