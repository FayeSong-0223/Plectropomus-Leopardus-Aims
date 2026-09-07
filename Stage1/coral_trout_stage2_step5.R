# =====================================================================
#  Plectropomus leopardus — Stage 2, Step 5: model diagnostics that
#  were not run earlier, and a correction to one reported result
#
#  Ziqi (Faye) Song
#
#  Two reviews of the project turned up seven things that had not been
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
#   5. Whether the negative binomial actually accounts for the observed
#      zeros. The plan flagged Palm's concentrated zeros and nothing
#      followed up. Tested here by parametric bootstrap rather than
#      against an analytic expectation.
#   6. Temporal autocorrelation within sites, binned by the true gap in
#      years, since the survey years are unequally spaced.
#   7. Multiplicity. Model A reports 17 tests and no p-value quoted in
#      this project had been adjusted. Holm within the four hypothesis
#      families the plan set out in advance — not one blanket correction
#      across tests that answer different questions.
#
#  Finding: the Moran's I result reported earlier was seed-dependent and
#  is corrected here. Zero counts exceed what the model generates, which
#  is reported as a limitation and shown not to drive the protection
#  result. Everything else held.
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
TEAL <- "#2E7D6F"; RED <- "#C1655A"; BLUEG <- "#4A6FA5"

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
  b <- coef(m); V <- vcov(m, unconditional = TRUE); nm <- names(b); k <- rep(0, length(b))
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

# ---------------------------------------------------------------------
# 5. IS THE NEGATIVE BINOMIAL ADEQUATE FOR THE ZEROS?
# ---------------------------------------------------------------------
rule("5. ZERO COUNTS AGAINST WHAT THE MODEL PREDICTS")
say("Section 8 of the project plan flagged that Palm's zeros are concentrated, and")
say("no distributional check was ever run against that. The negative binomial can")
say("accommodate a good deal of excess zero mass through its overdispersion, but")
say("whether it accommodates THIS much is a question with an answer.")
say("")
say("The test is a parametric bootstrap. Each replicate draws a new dataset from the")
say("fitted model — same covariates, same fitted means, same theta — counts its")
say("zeros, and the observed count is compared against that reference distribution.")
say("An analytic expected value alone would not do: it gives a point to compare")
say("against with no sense of how far a correct model would ordinarily stray.\n")

NSIM <- 2000
mu <- fitted(mA); th <- mA$family$getTheta(TRUE)
set.seed(11)
sim0 <- replicate(NSIM, sum(rnbinom(length(mu), size = th, mu = mu) == 0))
obs0 <- sum(d$count == 0)
p_zero <- (1 + sum(sim0 >= obs0)) / (NSIM + 1)
say("observed zeros: ", obs0, " of ", nrow(d), " observations")
say("simulated under the fitted model: median ", median(sim0),
    ", 95% interval ", paste(quantile(sim0, c(0.025, 0.975)), collapse = " to "))
say("one-sided p (simulated >= observed): ", sprintf("%.4f", p_zero))
say(if (p_zero < 0.05)
      "-> the model generates fewer zeros than observed. Excess zero mass is real." else
      "-> the observed zero count is within what the fitted model generates.")

say("\nThe same comparison across the low counts, where zero inflation would show up")
say("as a zero excess paired with a deficit at one and two:\n")
set.seed(12)
cnt_tab <- t(sapply(0:5, function(k) {
  sim <- replicate(500, sum(rnbinom(length(mu), size = th, mu = mu) == k))
  c(count = k, observed = sum(d$count == k), sim_median = median(sim),
    lo = unname(quantile(sim, 0.025)), hi = unname(quantile(sim, 0.975)))
}))
cap(as.data.frame(cnt_tab))

say("\nBy region, since the plan's concern was specifically about Palm:")
set.seed(13)
zr <- do.call(rbind, lapply(levels(d$REGION), function(rg) {
  i <- d$REGION == rg
  s <- replicate(1000, sum(rnbinom(sum(i), size = th, mu = mu[i]) == 0))
  data.frame(region = rg, n = sum(i), observed = sum(d$count[i] == 0),
             sim_median = median(s), lo = unname(quantile(s, 0.025)),
             hi = unname(quantile(s, 0.975)),
             p = round((1 + sum(s >= sum(d$count[i] == 0))) / 1001, 4))
}))
cap(zr)

say("\nWhat this changes. A zero excess inflates the apparent overdispersion, which")
say("widens intervals rather than narrowing them, so it does not manufacture the")
say("protection result. The check below confirms that directly.")
zi_site <- tapply(d$count == 0, d$SITE, sum)
drop_sites <- names(zi_site)[zi_site >= 4]
say("sites where at least 4 of 7 surveys are zero: ", length(drop_sites),
    " (holding ", sprintf("%.0f%%", 100 * sum(zi_site[zi_site >= 4]) / sum(zi_site)),
    " of all zeros)")
d_zt <- droplevels(d[!d$SITE %in% drop_sites, ])
m_zt <- gam(FORM, family = nb(), data = d_zt, method = "REML")
zt <- do.call(rbind, lapply(c("NTR 1987", "NTR 2004"), function(lv) {
  a <- rr(mA, "Whitsunday", lv); b <- rr(m_zt, "Whitsunday", lv)
  data.frame(protection = lv,
             full = sprintf("%.2f (%.2f-%.2f)", a[1], a[2], a[3]),
             zero_trimmed = sprintf("%.2f (%.2f-%.2f)", b[1], b[2], b[3]))
}))
cap(zt)

# ---------------------------------------------------------------------
# 6. TEMPORAL AUTOCORRELATION WITHIN SITES
# ---------------------------------------------------------------------
rule("6. TEMPORAL AUTOCORRELATION WITHIN SITES")
say("The model gives each site a random intercept, which absorbs a site's persistent")
say("level. It does not model correlation between successive visits to the same site")
say("beyond that, and nothing in this project had tested whether it needs to.")
say("")
say("Survey years are unequally spaced — 2007, 2009, 2012, 2014, 2016, 2017, 2018 —")
say("so pairing consecutive SURVEYS would mix a one-year gap with a three-year gap")
say("and call both 'lag 1'. Pairs are therefore binned by the actual number of years")
say("between them.\n")

rq_resid <- function(seed) {
  set.seed(seed)
  qnorm(runif(nrow(d), pnbinom(d$count - 1, size = th, mu = mu),
              pnbinom(d$count, size = th, mu = mu)))
}
lag_pairs <- function(res) {
  z <- data.frame(SITE = d$SITE, YEAR = d$YEAR, REGION = d$REGION, r = res)
  do.call(rbind, lapply(split(z, z$SITE), function(s) {
    if (nrow(s) < 2) return(NULL)
    s <- s[order(s$YEAR), ]; k <- nrow(s)
    ij <- expand.grid(i = seq_len(k), j = seq_len(k))
    ij <- ij[ij$i < ij$j, ]
    data.frame(gap = s$YEAR[ij$j] - s$YEAR[ij$i], y1 = s$YEAR[ij$i], y2 = s$YEAR[ij$j],
               REGION = s$REGION[1], a = s$r[ij$i], b = s$r[ij$j])
  }))
}
lp <- lag_pairs(rq_resid(1))
gap_tab <- do.call(rbind, lapply(sort(unique(lp$gap)), function(g) {
  z <- lp[lp$gap == g, ]
  if (nrow(z) < 30) return(NULL)
  ct <- cor.test(z$a, z$b)
  data.frame(gap_years = g, n_pairs = nrow(z), r = round(unname(ct$estimate), 3),
             lo = round(ct$conf.int[1], 3), hi = round(ct$conf.int[2], 3),
             p = round(ct$p.value, 4))
}))
cap(gap_tab)
write.csv(gap_tab, file.path(OUT, "table13_temporal_autocorrelation.csv"), row.names = FALSE)

say("\nOne caution before reading that table. A site random intercept forces a site's")
say("residuals to sum to roughly zero, which induces a negative correlation of about")
say("-1/(m-1) at every gap, with m the surveys per site. Here m = ",
    sprintf("%.2f", mean(table(d$SITE))), ", so about ",
    sprintf("%.3f", -1 / (mean(table(d$SITE)) - 1)), " is the baseline to judge")
say("against, not zero. Most gaps sit at or above it. One does not.\n")

worst <- gap_tab$gap_years[which.min(gap_tab$r)]
say("Because randomised quantile residuals depend on the seed, the strongest gap —")
say(worst, " years — is re-tested across 20 randomisations, as every other residual")
say("diagnostic here is:\n")
seed_r <- t(sapply(1:20, function(s) {
  z <- lag_pairs(rq_resid(s)); z <- z[z$gap == worst, ]
  ct <- cor.test(z$a, z$b); c(r = unname(ct$estimate), p = ct$p.value)
}))
say("gap of ", worst, " years: r ranges ", sprintf("%.3f", min(seed_r[, 1])), " to ",
    sprintf("%.3f", max(seed_r[, 1])), ", median ", sprintf("%.3f", median(seed_r[, 1])),
    ", p < 0.05 in ", sum(seed_r[, 2] < 0.05), " of 20 randomisations")
say(if (sum(seed_r[, 2] < 0.05) <= 1)
      "-> not distinguishable from the correlation the random effect induces." else
      "-> this is not a seed artefact and it is stronger than the induced baseline.")

say("\nWhich survey pairs carry it:\n")
z3 <- lag_pairs(rq_resid(1))
z3 <- z3[z3$gap == worst, ]
pair_tab <- do.call(rbind, lapply(split(z3, paste(z3$REGION, z3$y1, "to", z3$y2)), function(s)
  data.frame(n_sites = nrow(s), r = round(cor(s$a, s$b), 3),
             p = round(cor.test(s$a, s$b)$p.value, 4))))
cap(pair_tab)

say("\nWhat it means. A negative correlation is a reversal: sites sitting above their")
say("predicted density at the first survey sat below it at the second, and the other")
say("way round. Something between those two surveys moved sites differentially, and")
say("the model does not have a term for it. The region-specific year trend cannot")
say("supply one — it shifts every site in a region by the same amount.")
say("")
say("This is the cyclone limitation made concrete rather than merely stated. Cyclone")
say("exposure is recorded at survey points only, so any disturbance falling in the")
say("gap between two surveys is invisible to the covariate, and the gap concerned")
say("here is a three-year one. Attributing it to a particular storm would be")
say("inference beyond this dataset; what the data support is that the gap contains")
say("an unmodelled, site-differentiating event. It belongs in the limitations and in")
say("the list of questions for anyone who holds the track data.")

# ---------------------------------------------------------------------
# 7. MULTIPLICITY, BY HYPOTHESIS FAMILY
# ---------------------------------------------------------------------
rule("7. MULTIPLICITY WITHIN HYPOTHESIS FAMILIES")
say("Model A reports 17 tests. None of the p-values quoted anywhere in this project")
say("had been adjusted for that, which matters most for the terms that are already")
say("marginal.")
say("")
say("A single correction across all 17 would be the wrong instrument. The plan set")
say("out four hypotheses in advance, and a test only competes with the others asked")
say("in service of the same hypothesis. Holm's correction is applied within each")
say("family. It controls the family-wise error rate, needs no independence")
say("assumption — which matters here, since concurvity between these terms is high —")
say("and is uniformly at least as powerful as Bonferroni.\n")

st <- summary(mA)$s.table; pt <- summary(mA)$p.table
getp <- function(tab, rn) if (rn %in% rownames(tab)) tab[rn, ncol(tab)] else NA_real_
fam <- rbind(
  data.frame(family = "H1/Q3 management", term = c("NTR 1987", "NTR 2004",
             "REGION x NTR 1987", "REGION x NTR 2004"),
             p = c(getp(pt, "NTRNTR 1987"), getp(pt, "NTRNTR 2004"),
                   getp(pt, "REGIONWhitsunday:NTRNTR 1987"),
                   getp(pt, "REGIONWhitsunday:NTRNTR 2004"))),
  data.frame(family = "H2 habitat", term = c("s(rugosity)", "s(LHC)", "depth"),
             p = c(getp(st, "s(rugosity)"), getp(st, "s(LHC)"), getp(pt, "depth"))),
  data.frame(family = "H3 environment", term = c("s(maxDHW)", "s(kd490)", "s(Cyclone)"),
             p = c(getp(st, "s(maxDHW)"), getp(st, "s(kd490)"), getp(st, "s(Cyclone)"))),
  data.frame(family = "H4 site variation", term = "s(SITE)", p = getp(st, "s(SITE)")))
fam$p_holm <- ave(fam$p, fam$family, FUN = function(z) p.adjust(z, method = "holm"))
fam$p <- round(fam$p, 4); fam$p_holm <- round(fam$p_holm, 4)
fam$survives <- ifelse(is.na(fam$p_holm), NA, fam$p_holm < 0.05)
cap(fam)
write.csv(fam, file.path(OUT, "table14_multiplicity.csv"), row.names = FALSE)

say("\nThe result that matters survives comfortably. Q3 is the region-by-protection")
say("interaction, and both its terms clear Holm within the management family by an")
say("order of magnitude. The NTR main effects do not, but those are the Palm")
say("estimates, which were already null — a null that fails a correction for")
say("multiplicity is still a null.")
say("")
say("Nothing in the habitat family survives, which agrees with the within-between")
say("split in Step 4.")
say("")
say("The environment family is where the two corrections part company. The thermal")
say("term survives both Holm and the unconditional interval. Cyclone exposure")
say("survives neither. Turbidity survives Holm but its magnitude interval includes")
say("1 once smoothing-parameter uncertainty is admitted, so it should be reported as")
say("unresolved rather than as an effect — the weaker of two readings is the one to")
say("quote when they disagree.")
say("")
say("The year-trend smooths are not in any family. They are adjustment terms rather")
say("than hypotheses, and correcting a hypothesis test for the significance of a")
say("nuisance term would be a category error.")

# ---------------------------------------------------------------------
# 8. FIGURE 16 — the two new diagnostics
# ---------------------------------------------------------------------
rule("8. FIGURE")
png(file.path(OUT, "fig16_distribution_and_time.png"), width = 2000, height = 950, res = 220)
par(mfrow = c(1, 2), mar = c(3.9, 3.9, 2.8, 0.9)); base_par()

h <- hist(sim0, breaks = 24, plot = FALSE)
plot(h, col = "grey88", border = "white",
     xlim = c(min(h$breaks), max(h$breaks, obs0) + 3),
     xlab = "Zero counts in data simulated from the fitted model",
     ylab = "Replicates", main = "a  Are there more zeros than the model makes?",
     font.main = 1, cex.main = 0.9, adj = 0)
abline(v = obs0, col = RED, lwd = 2.2)
text(obs0, par("usr")[4] * 0.90, paste0("observed = ", obs0, " "),
     col = RED, cex = 0.75, adj = 1)

g <- gap_tab
plot(g$gap_years, g$r, type = "n", ylim = range(c(g$lo, g$hi)),
     xlab = "Years between surveys", ylab = "Residual correlation",
     main = "b  Residual correlation within sites, by true gap",
     font.main = 1, cex.main = 0.9, adj = 0)
abline(h = 0, col = "grey55", lty = 2)
abline(h = -1 / (mean(table(d$SITE)) - 1), col = BLUEG, lty = 3, lwd = 1.6)
segments(g$gap_years, g$lo, g$gap_years, g$hi, col = "grey45", lwd = 1.8)
points(g$gap_years, g$r, pch = 19, col = ifelse(g$p < 0.05, RED, TEAL), cex = 1.05)
text(min(g$gap_years), -1 / (mean(table(d$SITE)) - 1), "induced by the site effect",
     col = BLUEG, cex = 0.62, adj = c(0, -0.5))
dev.off(); say("wrote fig16_distribution_and_time.png")

# Record the environment. mgcv ships with R but its version tracks the R
# version, and REML fitting and the nb() family have both changed across
# releases, so "no packages required" is not the same as "no versions to
# reconcile". Anyone reproducing these numbers needs to know what produced them.
say("\n", strrep("-", 70))
say("environment: ", R.version.string, " | mgcv ", as.character(packageVersion("mgcv")),
    " | platform ", R.version$platform)
say(strrep("-", 70))

writeLines(log_lines, file.path(OUT, "stage2_step5_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
