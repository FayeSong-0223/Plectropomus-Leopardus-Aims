# =====================================================================
#  Plectropomus leopardus — Stage 2, Step 4: completing the plan
#  Palm and Whitsunday inshore reefs, Great Barrier Reef, 2007-2018
#
#  Ziqi (Faye) Song
#
#  An audit against the project plan found four specified deliverables
#  outstanding. This script produces them.
#
#   1. Effect magnitudes across the observed range of each covariate.
#      Q1 asks "how large is each association across the observed range
#      of that variable" — a question p-values and edf do not answer.
#   2. Moran's I on site-level residuals (plan section 10, check 4).
#      Tests whether the site random effect has absorbed the spatial
#      structure, using the coordinates parsed in Stage 1.
#   3. The variable dictionary (plan section 11, Table 1) — definition,
#      scale of variation, retained or dropped, and why.
#   4. The shared components of the variance decomposition, shown in the
#      figure rather than only in the log (plan section 11, Fig 4).
#   5. The within-between (Mundlak) split on the habitat block, which the
#      plan's section 9 scopes specifically to Q4. A pooled habitat
#      coefficient conflates 'sites with more structure hold more fish'
#      with 'a site that gains structure gains fish'. Q4 asks the second.
#
#  Dependencies: mgcv only.
#  Run:  Rscript coral_trout_stage2_step4.R
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
TEAL <- "#2E7D6F"; RED <- "#C1655A"; BLUE <- "#4A6FA5"

# ---------------------------------------------------------------------
# 0. REFIT MODEL A
# ---------------------------------------------------------------------
rule("0. MODEL")
f <- file.path(OUT, "analysis_dataset.csv")
if (!file.exists(f)) stop("outputs/analysis_dataset.csv not found. Run coral_trout_stage1.R first.")
d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
names(d)[names(d) == "LHC_%"] <- "LHC"; names(d)[names(d) == "LCC_%"] <- "LCC"
names(d)[names(d) == "Corrected depth"] <- "depth"
d$REGION   <- factor(d$REGION,   levels = c("Palm", "Whitsunday"))
d$NTR      <- factor(d$NTR,      levels = c("Fished", "NTR 1987", "NTR 2004"))
d$EXPOSURE <- factor(d$EXPOSURE, levels = c("Sheltered", "Semi-Exposed", "Exposed"))
d$SITE     <- factor(d$SITE); d$RY <- factor(paste(d$REGION, d$YEAR))

mA <- gam(count ~ REGION * NTR + EXPOSURE + depth + s(YEAR, by = REGION, k = 5) +
            s(rugosity, k = 5) + s(LHC, k = 5) +
            s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) +
            s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")
say("Model A refitted: deviance explained ", sprintf("%.1f%%", summary(mA)$dev.expl * 100))

# ---------------------------------------------------------------------
# 1. EFFECT MAGNITUDE ACROSS THE OBSERVED RANGE
# ---------------------------------------------------------------------
rule("1. HOW LARGE IS EACH ASSOCIATION?")
say("Q1 asks how large each association is across the observed range of the")
say("variable. A p-value does not answer that and neither does an edf. Below,")
say("each covariate is moved from its 10th to its 90th percentile with everything")
say("else held fixed, and the change in predicted density is reported as a ratio.")
say("Because the contrast is a difference of two rows of the same design matrix,")
say("every term that does not involve the covariate cancels exactly — including")
say("the site random effect.\n")

eff_range <- function(m, v, dat, lo = 0.10, hi = 0.90) {
  q <- quantile(dat[[v]], c(lo, hi), na.rm = TRUE)
  nd <- dat[c(1, 1), ]; nd[[v]] <- as.numeric(q)
  X  <- predict(m, newdata = nd, type = "lpmatrix")
  dX <- X[2, ] - X[1, ]
  est <- sum(dX * coef(m))
  # Unconditional: includes uncertainty in the smoothing parameters. Three of the
  # four intervals that excluded 1 under the conditional matrix no longer do.
  se  <- sqrt(as.numeric(t(dX) %*% vcov(m, unconditional = TRUE) %*% dX))
  c(from = q[[1]], to = q[[2]], ratio = exp(est),
    lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

covs <- c("rugosity", "LHC", "depth", "kd490", "maxDHW", "Cyclone")
nice <- c("Rugosity index", "Live hard coral (%)", "Depth (m)",
          "kd490 (turbidity)", "Max degree heating weeks", "Cyclone exposure index")
eff <- do.call(rbind, lapply(seq_along(covs), function(i) {
  v <- eff_range(mA, covs[i], d)
  data.frame(covariate = nice[i], p10 = round(v[["from"]], 3), p90 = round(v[["to"]], 3),
             ratio = round(v[["ratio"]], 3), lo = round(v[["lo"]], 3), hi = round(v[["hi"]], 3),
             pct_change = sprintf("%+.0f%%", (v[["ratio"]] - 1) * 100))
}))
cap(eff)
write.csv(eff, file.path(OUT, "table8_effect_magnitudes.csv"), row.names = FALSE)
say("\nIntervals containing 1 mean the association is not distinguishable from no")
say("effect across the range the data actually cover.")

# protection, for comparison, on the same scale
rr <- function(m, region, level) {
  b <- coef(m); V <- vcov(m, unconditional = TRUE); nm <- names(b); k <- rep(0, length(b))
  k[match(paste0("NTR", level), nm)] <- 1
  ix <- paste0("REGIONWhitsunday:NTR", level)
  if (region == "Whitsunday") k[match(ix, nm)] <- 1
  e <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(exp(e), exp(e - 1.96 * se), exp(e + 1.96 * se))
}
prot <- do.call(rbind, lapply(c("Palm", "Whitsunday"), function(r)
  do.call(rbind, lapply(c("NTR 1987", "NTR 2004"), function(l) {
    v <- rr(mA, r, l)
    data.frame(covariate = paste0(r, ", ", l, " vs fished"),
               p10 = NA, p90 = NA, ratio = round(v[1], 3),
               lo = round(v[2], 3), hi = round(v[3], 3),
               pct_change = sprintf("%+.0f%%", (v[1] - 1) * 100))
  }))))
say("\nprotection contrasts on the same ratio scale, for comparison:")
cap(prot)

# ---------------------------------------------------------------------
# 2. SPATIAL AUTOCORRELATION IN RESIDUALS  (plan section 10, check 4)
# ---------------------------------------------------------------------
rule("2. SPATIAL RESIDUAL CHECK")
say("If the site random effect is doing its job, site-level mean residuals should")
say("show no spatial structure. Sites within an island group sit a few hundred")
say("metres apart, so this is tested WITHIN each region: across regions the only")
say("thing a global test would detect is the 200 km gap between them.\n")

qres <- function(m, y) {
  mu <- fitted(m); th <- m$family$getTheta(TRUE)
  qnorm(runif(length(y), pnbinom(y - 1, size = th, mu = mu), pnbinom(y, size = th, mu = mu)))
}
d$r <- qres(mA, d$count)
site <- aggregate(cbind(r, lat_dd, long_dd) ~ SITE + REGION, data = d, FUN = mean)

moran <- function(x, lon, lat, nperm = 4999) {
  n <- length(x); R <- 6371; latm <- mean(lat) * pi / 180
  ex <- R * (lon * pi / 180) * cos(latm); ey <- R * (lat * pi / 180)
  D <- as.matrix(dist(cbind(ex, ey)))
  W <- 1 / D; diag(W) <- 0; W[!is.finite(W)] <- 0
  z <- x - mean(x); S <- sum(W)
  I <- (n / S) * sum(W * outer(z, z)) / sum(z^2)
  perm <- replicate(nperm, { zz <- sample(z); (n / S) * sum(W * outer(zz, zz)) / sum(zz^2) })
  c(I = I, expected = -1 / (n - 1), p = (sum(abs(perm) >= abs(I)) + 1) / (nperm + 1),
    med_km = median(D[upper.tri(D)]), min_km = min(D[upper.tri(D)]))
}
mor <- do.call(rbind, lapply(levels(site$REGION), function(rg) {
  s <- site[site$REGION == rg, ]
  v <- moran(s$r, s$long_dd, s$lat_dd)
  data.frame(region = rg, n_sites = nrow(s), Moran_I = round(v[["I"]], 4),
             expected = round(v[["expected"]], 4), p_perm = round(v[["p"]], 4),
             median_sep_km = round(v[["med_km"]], 2), min_sep_km = round(v[["min_km"]], 3))
}))
cap(mor)
say("\nInverse-distance weights, 4999 permutations. Under no spatial structure")
say("Moran's I sits near -1/(n-1), which is slightly below zero, not at zero.")
say("\nCAUTION — these are SINGLE-SEED values and must not be read as a result.")
say("Randomised quantile residuals draw from a uniform, so Moran's I and its")
say("permutation p-value both depend on set.seed(). Step 5 repeats this across")
say("twenty seeds and supersedes the numbers above. Its finding: Whitsunday's I is")
say("positive under every seed (0.079 to 0.167) while Palm's is near zero under")
say("every seed, but the Whitsunday permutation p falls below 0.05 in only 11 of")
say("20 randomisations. Residual spatial structure in Whitsunday is suggested,")
say("not established. See coral_trout_stage2_step5.R and stage2_step5_log.txt.")
write.csv(mor, file.path(OUT, "table9_spatial_autocorrelation.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 3. VARIABLE DICTIONARY  (plan section 11, Table 1)
# ---------------------------------------------------------------------
rule("3. VARIABLE DICTIONARY")
ry <- paste(d$REGION, d$YEAR)
bshare <- function(v) { if (is.null(d[[v]])) return(NA); x <- d[[v]]; m <- ave(x, ry, FUN = function(z) mean(z, na.rm = TRUE))
                        round(var(m, na.rm = TRUE) / var(x, na.rm = TRUE), 3) }
wsite  <- function(v) if (is.null(d[[v]])) NA else
  sum(tapply(d[[v]], d$SITE, function(z) length(unique(z))) > 1)
dict <- rbind(
 data.frame(variable="pms.leop / count", role="Response", varies_at="Site-year",
   between_RY_share=bshare("pms.leop"), status="Retained",
   reason="Density recovered to integer counts; area constant so no offset"),
 data.frame(variable="NTR", role="Management", varies_at="Site (fixed)", between_RY_share=NA,
   status="Retained", reason="Three levels kept separate; NTR Pooled not used"),
 data.frame(variable="EXPOSURE", role="Confounder", varies_at="Site (fixed)", between_RY_share=NA,
   status="Retained", reason="Confounded with protection; must accompany it"),
 data.frame(variable="REGION, YEAR", role="Structure", varies_at="-", between_RY_share=NA,
   status="Retained", reason="Year made region-specific after Step 1 residual diagnostics"),
 data.frame(variable="SITE", role="Structure", varies_at="-", between_RY_share=NA,
   status="Retained", reason="Random intercept, 71 levels, REML"),
 data.frame(variable="rugosity", role="Habitat", varies_at="Site-year",
   between_RY_share=bshare("rugosity"), status="Retained",
   reason="Structural complexity; chosen over SCI (r = 0.93)"),
 data.frame(variable="LHC_%", role="Habitat", varies_at="Site-year",
   between_RY_share=bshare("LHC"), status="Retained",
   reason="Live hard coral; chosen over LCC_% (r = 0.80) in advance"),
 data.frame(variable="Corrected depth", role="Habitat", varies_at=paste0("Varies at ", wsite("depth"), " of 71 sites"),
   between_RY_share=bshare("depth"), status="Retained", reason="Near site-constant in practice"),
 data.frame(variable="kd490", role="Environment", varies_at="Site-year",
   between_RY_share=bshare("kd490"), status="Retained",
   reason="Water clarity; chosen over ChlA (r = 1.00 in this subset)"),
 data.frame(variable="maxDHW", role="Disturbance", varies_at="Region-year",
   between_RY_share=bshare("maxDHW"), status="Retained (Model A only)",
   reason="Nearly nested in region-year; not identifiable in Model B"),
 data.frame(variable="Cyclone", role="Disturbance", varies_at="Mostly site-level",
   between_RY_share=bshare("Cyclone"), status="Retained",
   reason="Better identified than maxDHW; interval integration impossible from this extract"),
 data.frame(variable="SCI", role="Habitat", varies_at="Site-year",
   between_RY_share=bshare("SCI"), status="Dropped", reason="Duplicate of rugosity (r = 0.93)"),
 data.frame(variable="ChlA", role="Environment", varies_at="Site-year",
   between_RY_share=bshare("ChlA"), status="Dropped", reason="Duplicate of kd490 (r = 1.00)"),
 data.frame(variable="LCC_%", role="Habitat", varies_at="Site-year",
   between_RY_share=bshare("LCC"), status="Dropped",
   reason="Correlates 0.80 with LHC_%; NOT reinstated after the habitat null, which would be selection on the outcome"),
 data.frame(variable="SSTmean", role="Environment", varies_at="Region-year",
   between_RY_share=bshare("SSTmean"), status="Dropped",
   reason="98% between-region-year; would compete with maxDHW for the same dozen contrasts"),
 data.frame(variable="LT Fprimary", role="Fishing proxy", varies_at="Site-constant (0 of 71 sites vary)",
   between_RY_share=NA, status="Dropped",
   reason="Modelled index not a measurement, and site-constant, so competes with protection and exposure for the same 71 between-site df"),
 data.frame(variable="wave exposure index", role="Confounder", varies_at="Site-constant (0 of 71 sites vary)",
   between_RY_share=NA, status="Dropped",
   reason="Redundant with categorical EXPOSURE; estimable but not needed alongside it"))
cap(dict)
write.csv(dict, file.path(OUT, "table1_variable_dictionary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4. VARIANCE DECOMPOSITION WITH SHARED COMPONENTS SHOWN
# ---------------------------------------------------------------------
rule("4. VARIANCE DECOMPOSITION, UNIQUE AND SHARED")
tm <- predict(mA, type = "terms"); cn <- colnames(tm)
blk <- list(Management = grep("NTR", cn, value = TRUE),
            Habitat    = unique(c(grep("s\\((rugosity|LHC)\\)", cn, value = TRUE), grep("^depth$", cn, value = TRUE))),
            Environment= grep("s\\((kd490|maxDHW|Cyclone)\\)", cn, value = TRUE),
            `Space/time`= unique(c(grep("^REGION$|^EXPOSURE$", cn, value = TRUE),
                                   grep("YEAR", cn, value = TRUE))),
            `Site (RE)`= grep("s\\(SITE\\)", cn, value = TRUE))
blk <- lapply(blk, function(x) unique(x[x %in% cn]))
B <- sapply(blk, function(k) if (length(k)) rowSums(tm[, k, drop = FALSE]) else rep(0, nrow(tm)))
S <- cov(B); tot <- sum(S); uni <- diag(S) / tot
sh <- S; diag(sh) <- NA
shd <- as.data.frame(as.table(round(2 * sh / tot, 4)))
shd <- shd[!is.na(shd$Freq) & as.character(shd$Var1) < as.character(shd$Var2), ]
shd <- shd[order(-abs(shd$Freq)), ]; names(shd) <- c("block_1", "block_2", "shared_share")
say("unique shares:"); cap(round(sort(uni, decreasing = TRUE), 3))
say("\nshared shares (2 x covariance / total, so unique + shared sums to 1):")
cap(shd)
say("\ncheck: unique ", sprintf("%.3f", sum(uni)), " + shared ",
    sprintf("%.3f", sum(shd$shared_share)), " = ",
    sprintf("%.3f", sum(uni) + sum(shd$shared_share)))

png(file.path(OUT, "fig7_variance_decomposition.png"), width = 2100, height = 1150, res = 220)
par(mfrow = c(1, 2), mar = c(3.9, 7.6, 2.6, 1.2)); base_par()
o <- sort(uni)
bp <- barplot(o, horiz = TRUE, col = TEAL, border = NA, xlim = c(0, max(o) * 1.3),
              xlab = "Share of variance", main = "a  Unique contribution",
              font.main = 1, cex.main = 0.98, adj = 0, las = 1, cex.names = 0.78)
text(o, bp, sprintf(" %.3f", o), pos = 4, cex = 0.72, col = "grey25", xpd = NA)
s2 <- shd[order(shd$shared_share), ]
lab2 <- paste(s2$block_1, "+", s2$block_2)
bp2 <- barplot(s2$shared_share, horiz = TRUE, names.arg = lab2,
               col = ifelse(s2$shared_share > 0, BLUE, RED), border = NA,
               xlim = range(c(s2$shared_share, 0)) * 1.45,
               xlab = "Share of variance", main = "b  Shared between blocks",
               font.main = 1, cex.main = 0.98, adj = 0, las = 1, cex.names = 0.62)
abline(v = 0, col = "grey45");
text(s2$shared_share, bp2, sprintf("%.3f", s2$shared_share),
     pos = ifelse(s2$shared_share > 0, 4, 2), cex = 0.66, col = "grey25", xpd = NA)
dev.off(); say("\nwrote fig7_variance_decomposition.png (now two panels)")

## Figure 12 — effect magnitudes
png(file.path(OUT, "fig12_effect_magnitudes.png"), width = 2000, height = 1150, res = 220)
par(mar = c(3.9, 9.8, 2.6, 1.2)); base_par()
ae <- rbind(eff[, c("covariate", "ratio", "lo", "hi")], prot[, c("covariate", "ratio", "lo", "hi")])
ord <- rev(seq_len(nrow(ae)))
plot(NA, xlim = range(c(ae$lo, ae$hi, 1)), ylim = c(0.5, nrow(ae) + 0.5), log = "x", yaxt = "n",
     xlab = "Density ratio (10th to 90th percentile, or vs fished)", ylab = "",
     main = "How large is each association across the range the data cover?",
     font.main = 1, cex.main = 0.92, adj = 0)
axis(2, at = ord, labels = ae$covariate, cex.axis = 0.72)
abline(v = 1, col = "grey60", lty = 2)
col <- c(rep(TEAL, nrow(eff)), rep(BLUE, nrow(prot)))
segments(ae$lo, ord, ae$hi, ord, col = col, lwd = 2.2)
points(ae$ratio, ord, pch = 19, col = col, cex = 1.1)
legend("bottomleft", legend = c("covariate, 10th to 90th percentile", "protection vs fished"),
       col = c(TEAL, BLUE), lwd = 2.2, pch = 19, bty = "n", cex = 0.72, text.col = "grey20")
dev.off(); say("wrote fig12_effect_magnitudes.png")

## Figure 13 — spatial residual check
png(file.path(OUT, "fig13_spatial_residuals.png"), width = 2000, height = 1000, res = 220)
par(mfrow = c(1, 2), mar = c(3.8, 3.9, 2.6, 0.9)); base_par()
for (rg in levels(site$REGION)) {
  s <- site[site$REGION == rg, ]
  R <- 6371; latm <- mean(s$lat_dd) * pi / 180
  ex <- R * (s$long_dd * pi / 180) * cos(latm); ey <- R * (s$lat_dd * pi / 180)
  D <- as.matrix(dist(cbind(ex, ey))); z <- s$r - mean(s$r)
  prod <- outer(z, z); ut <- upper.tri(D)
  mi <- mor[mor$region == rg, ]
  plot(D[ut], prod[ut], pch = 16, col = adjustcolor("grey30", 0.35), cex = 0.55,
       xlab = "Distance between sites (km)", ylab = "Product of residual deviations",
       main = sprintf("%s   Moran's I = %.3f  (single seed - see Fig 15)", rg, mi$Moran_I),
       font.main = 1, cex.main = 0.85, adj = 0)
  abline(h = 0, col = "grey55", lty = 2)
  lines(lowess(D[ut], prod[ut]), col = RED, lwd = 2)
}
dev.off(); say("wrote fig13_spatial_residuals.png")

# ---------------------------------------------------------------------
# 5. WITHIN-BETWEEN (MUNDLAK) SPLIT ON THE HABITAT BLOCK
# ---------------------------------------------------------------------
rule("5. WITHIN AND BETWEEN HABITAT EFFECTS (Q4)")
say("A pooled habitat coefficient answers two questions at once and separates")
say("neither: do sites with more structure hold more fish (between-site, and")
say("confounded by everything else fixed about a site), and does a site that gains")
say("structure gain fish (within-site, immune to any time-invariant site")
say("characteristic). Q4 asks the second. Splitting each habitat covariate into a")
say("site mean and a within-site deviation estimates both.\n")

for (v in c("rugosity", "LHC")) {
  d[[paste0(v, "_bw")]] <- ave(d[[v]], d$SITE, FUN = function(z) mean(z, na.rm = TRUE))
  d[[paste0(v, "_wi")]] <- d[[v]] - d[[paste0(v, "_bw")]]
}
say("depth is site-constant in this subset (varies at 0 of 71 sites), so it has no")
say("within-site component and is left unsplit.\n")

mM <- gam(count ~ REGION * NTR + EXPOSURE + depth + s(YEAR, by = REGION, k = 5) +
            s(rugosity_bw, k = 5) + s(rugosity_wi, k = 5) +
            s(LHC_bw, k = 5) + s(LHC_wi, k = 5) +
            s(kd490, k = 5) + s(maxDHW, k = 5) +
            s(Cyclone, k = 5) + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")
say("Mundlak model: deviance explained ", sprintf("%.1f%%", summary(mM)$dev.expl * 100),
    "  (pooled model: ", sprintf("%.1f%%", summary(mA)$dev.expl * 100), ")")
say("\nsmooth terms:")
cap(round(summary(mM)$s.table[grep("rugosity|LHC", rownames(summary(mM)$s.table)), , drop = FALSE], 4))

mund <- do.call(rbind, lapply(
  list(c("rugosity_bw", "Rugosity, between sites"), c("rugosity_wi", "Rugosity, within site"),
       c("LHC_bw", "Live hard coral, between sites"), c("LHC_wi", "Live hard coral, within site")),
  function(z) {
    v <- eff_range(mM, z[1], d)
    data.frame(term = z[2], p10 = round(v[["from"]], 3), p90 = round(v[["to"]], 3),
               ratio = round(v[["ratio"]], 3), lo = round(v[["lo"]], 3), hi = round(v[["hi"]], 3))
  }))
say("\nmagnitude across the 10th-90th percentile of each component:")
cap(mund)
write.csv(mund, file.path(OUT, "table10_within_between_habitat.csv"), row.names = FALSE)

rug_wi <- mund[mund$term == "Rugosity, within site", ]
rug_bw <- mund[mund$term == "Rugosity, between sites", ]
say("\nReading for Q4. The within-site rugosity estimate is the one Q4 asks for,")
say("because it is immune to confounding by any fixed characteristic of a site.")
say("Within-site: ", sprintf("%.2f (%.2f-%.2f)", rug_wi$ratio, rug_wi$lo, rug_wi$hi),
    "   Between-site: ", sprintf("%.2f (%.2f-%.2f)", rug_bw$ratio, rug_bw$lo, rug_bw$hi))
say("If these disagree, the pooled estimate reported earlier was averaging two")
say("different things and should not be quoted without the split.")

png(file.path(OUT, "fig14_within_between_habitat.png"), width = 1900, height = 1000, res = 220)
par(mar = c(3.9, 10.6, 2.6, 1.2)); base_par()
ord <- rev(seq_len(nrow(mund)))
plot(NA, xlim = range(c(mund$lo, mund$hi, 1)), ylim = c(0.5, nrow(mund) + 0.5), log = "x",
     yaxt = "n", xlab = "Density ratio, 10th to 90th percentile", ylab = "",
     main = "Habitat effects split into between-site and within-site components",
     font.main = 1, cex.main = 0.92, adj = 0)
axis(2, at = ord, labels = mund$term, cex.axis = 0.74)
abline(v = 1, col = "grey60", lty = 2)
col <- ifelse(grepl("within", mund$term), TEAL, "grey55")
segments(mund$lo, ord, mund$hi, ord, col = col, lwd = 2.2)
points(mund$ratio, ord, pch = 19, col = col, cex = 1.1)
legend("bottomleft", legend = c("within site (what Q4 asks)", "between sites"),
       col = c(TEAL, "grey55"), lwd = 2.2, pch = 19, bty = "n", cex = 0.74, text.col = "grey20")
dev.off(); say("\nwrote fig14_within_between_habitat.png")

# Record the environment. mgcv ships with R but its version tracks the R
# version, and REML fitting and the nb() family have both changed across
# releases, so "no packages required" is not the same as "no versions to
# reconcile". Anyone reproducing these numbers needs to know what produced them.
say("\n", strrep("-", 70))
say("environment: ", R.version.string, " | mgcv ", as.character(packageVersion("mgcv")),
    " | platform ", R.version$platform)
say(strrep("-", 70))

writeLines(log_lines, file.path(OUT, "stage2_step4_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
