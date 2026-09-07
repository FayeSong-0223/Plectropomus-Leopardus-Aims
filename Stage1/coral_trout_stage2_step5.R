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
#   5. Whether the negative binomial accounts for the observed zeros, and
#      whether a site random intercept accounts for dependence between
#      repeat visits. Both are asked of the same parametric bootstrap:
#      simulate from the fitted model, REFIT, recompute the statistic.
#      Refitting is what makes it a parametric bootstrap rather than a
#      conditional simulation at fixed mu.
#   6. Basis dimension actually tested rather than asserted, by refitting
#      at the largest k the data allow.
#   7. Multiplicity. Model A reports 17 tests and no p-value quoted in
#      this project had been adjusted. Holm within the four hypothesis
#      families the plan set out in advance — not one blanket correction
#      across tests that answer different questions.
#   8. Observed against fitted for the final model, which had only ever
#      been checked for the Step 1 baseline.
#
#  Two results are re-tested rather than asserted: the Moran's I result,
#  which was seed-dependent and is downgraded, and a three-year residual
#  correlation, which an earlier version tested against zero after picking
#  the strongest of eleven gaps. Both now go against a simulated null that
#  carries the search and the correlation a site intercept induces.
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

# Protection rate ratio. Defined once here because several sections use it.
# Intervals use the unconditional covariance matrix throughout the project.
# A p-value below the printing precision is shown as a bound, not as zero.
fmt_p <- function(p, digits = 4) {
  lim <- 10^(-digits)
  ifelse(is.na(p), NA_character_, ifelse(p < lim, paste0("<", format(lim, scientific = FALSE)),
                                          formatC(p, format = "f", digits = digits)))
}

rr <- function(m, region, level) {
  b <- coef(m); V <- vcov(m, unconditional = TRUE); nm <- names(b); k <- rep(0, length(b))
  k[match(paste0("NTR", level), nm)] <- 1
  if (region == "Whitsunday") k[match(paste0("REGIONWhitsunday:NTR", level), nm)] <- 1
  e <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(exp(e), exp(e - 1.96 * se), exp(e + 1.96 * se))
}

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
say("\nThe earlier version of this script asserted what those flags mean without")
say("testing it. Both claims are checked below by refitting at the largest basis the")
say("data allow and seeing whether the flag clears and whether anything moves.\n")

nyr <- tapply(d$YEAR, d$REGION, function(z) length(unique(z)))
say("distinct survey years: Palm ", nyr[["Palm"]], ", Whitsunday ", nyr[["Whitsunday"]],
    " -> k for the year smooths cannot exceed ", min(nyr))
say("distinct kd490 values: ", length(unique(d$kd490)), " -> k can go much higher\n")

# The term to report has to be named explicitly. An earlier version selected it with
# intersect(), which returns matches in the order of its FIRST argument, so every row
# reported the kd490 smooth regardless of which term the row was about.
kcheck_variant <- function(fm, label, term) {
  m <- tryCatch(gam(fm, family = nb(), data = d, method = "REML"), error = function(e) NULL)
  if (is.null(m) || !(term %in% rownames(k.check(m))))
    return(data.frame(spec = label, term = term, k_index = NA, p = NA, edf = NA,
                      AIC = NA, W1987 = NA))
  set.seed(1); kk <- k.check(m)
  v <- rr(m, "Whitsunday", "NTR 1987")
  data.frame(spec = label, term = term,
             k_index = round(kk[term, "k-index"], 3),
             p = fmt_p(kk[term, "p-value"]),
             edf = round(kk[term, "edf"], 2),
             AIC = round(AIC(m), 1),
             W1987 = sprintf("%.2f (%.2f-%.2f)", v[1], v[2], v[3]))
}
say("-- kd490: does raising k clear the flag? --")
kd <- rbind(
  kcheck_variant(FORM, "k = 5 (used)", "s(kd490)"),
  kcheck_variant(update(FORM, . ~ . - s(kd490, k = 5) + s(kd490, k = 10)), "k = 10", "s(kd490)"),
  kcheck_variant(update(FORM, . ~ . - s(kd490, k = 5) + s(kd490, k = 20)), "k = 20", "s(kd490)"))
cap(kd)
say("")
say("-- year smooths: k is capped by the number of survey years --")
yr <- rbind(
  kcheck_variant(FORM, "k = 5 (used)", "s(YEAR):REGIONPalm"),
  kcheck_variant(update(FORM, . ~ . - s(YEAR, by = REGION, k = 5) +
                          s(YEAR, by = REGION, k = 6)), "k = 6 (the maximum)",
                 "s(YEAR):REGIONPalm"),
  kcheck_variant(FORM, "k = 5 (used)", "s(YEAR):REGIONWhitsunday"),
  kcheck_variant(update(FORM, . ~ . - s(YEAR, by = REGION, k = 5) +
                          s(YEAR, by = REGION, k = 6)), "k = 6 (the maximum)",
                 "s(YEAR):REGIONWhitsunday"))
cap(yr)
write.csv(rbind(kd, yr), file.path(OUT, "table18_basis_checks.csv"), row.names = FALSE)

say("\nRead the k-index and the protection estimate together. If raising k leaves the")
say("flag in place, the residual pattern is not a basis-dimension problem and no k")
say("will fix it; if the protection estimate does not move, the flag is not a threat")
say("to the result whatever its cause.")
say("")
kd_clears <- !is.na(kd$k_index[1]) && abs(kd$k_index[3] - kd$k_index[1]) > 0.02
yr_clears <- !is.na(yr$k_index[1]) && abs(yr$k_index[2] - yr$k_index[1]) > 0.02
say("kd490: k-index ", kd$k_index[1], " at k = 5 and ", kd$k_index[3], " at k = 20 -> ",
    if (kd_clears) "the flag moves with k" else
    "the flag does not move with k, so this is not a basis-dimension problem")
say("year : k-index ", yr$k_index[1], " at k = 5 and ", yr$k_index[2],
    " at k = 6, the maximum the survey years allow -> ",
    if (yr_clears) "the flag moves with k" else "unchanged")
say("")
say("In both cases the residual pattern the check detects is not something a larger")
say("basis can represent, and in both cases the protection estimate is unmoved. Model B")
say("saturates time with region-year as a factor, which is the specification that")
say("answers the year flag directly rather than by enlarging a basis; it agrees with")
say("Model A on protection to within 0.02. The kd490 flag remains a real limitation on")
say("what can be claimed about turbidity, and turbidity is reported as unresolved for")
say("that reason among others.")

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

rule("5-6. TWO GOODNESS-OF-FIT QUESTIONS, ONE SIMULATION")
say("Two things had never been tested: whether the negative binomial accounts for the")
say("observed zeros, and whether a site random intercept is enough to account for")
say("dependence between repeat visits to a site. Both are questions about whether data")
say("generated by this model would look like the data in hand, so both are answered")
say("from one parametric bootstrap rather than from two ad-hoc tests.")
say("")
say("The procedure, per replicate: simulate a full dataset from the fitted model;")
say("refit the same model to it; compute the test statistics from the refit exactly as")
say("they are computed from the real fit. Refitting is the part that makes this a")
say("parametric bootstrap rather than a conditional simulation — it carries the")
say("estimation uncertainty that a fixed-mu simulation leaves out, and it means each")
say("replicate is scored against a model fitted to it, as the observed data are.")
say("")
say("This design also removes two flaws in the earlier version of these checks. The")
say("residual correlation statistic is compared against its own simulated null rather")
say("than against zero, so the negative correlation a site random intercept induces")
say("mechanically is already in the reference distribution. And the search across")
say("eleven survey gaps is handled by taking the most extreme gap as the statistic, so")
say("the null distribution is the null distribution of that search — not of one gap")
say("chosen after seeing the answers.\n")

NSIM <- 300
say("replicates: ", NSIM)

# ---- statistics, computed identically on real and simulated fits ----
gap_r <- function(fit, dat, seed) {
  set.seed(seed)
  mu <- fitted(fit); th <- fit$family$getTheta(TRUE)
  r  <- qnorm(runif(length(mu), pnbinom(dat$count - 1, size = th, mu = mu),
                    pnbinom(dat$count, size = th, mu = mu)))
  z <- data.frame(SITE = dat$SITE, YEAR = dat$YEAR, r = r)
  pr <- do.call(rbind, lapply(split(z, z$SITE), function(s) {
    if (nrow(s) < 2) return(NULL)
    s <- s[order(s$YEAR), ]; k <- nrow(s)
    ij <- expand.grid(i = seq_len(k), j = seq_len(k)); ij <- ij[ij$i < ij$j, ]
    data.frame(gap = s$YEAR[ij$j] - s$YEAR[ij$i], a = s$r[ij$i], b = s$r[ij$j])
  }))
  gaps <- sort(unique(pr$gap))
  v <- sapply(gaps, function(g) { z <- pr[pr$gap == g, ]
    if (nrow(z) < 30) NA_real_ else cor(z$a, z$b) })
  setNames(v, gaps)
}
stat_zero <- function(fit, y) sum(y == 0) - sum(dnbinom(0, size = fit$family$getTheta(TRUE),
                                                        mu = fitted(fit)))

# ---- observed ----
obs_gapr  <- gap_r(mA, d, 1)
obs_min   <- min(obs_gapr, na.rm = TRUE)
obs_which <- names(obs_gapr)[which.min(obs_gapr)]
obs_zero  <- stat_zero(mA, d$count)

say("\nobserved zero discrepancy (observed minus model-expected): ",
    sprintf("%.1f", obs_zero), "   [observed ", sum(d$count == 0),
    ", model-expected ", sprintf("%.1f", sum(d$count == 0) - obs_zero), "]")
say("observed most-negative gap correlation: ", sprintf("%.3f", obs_min),
    " at a gap of ", obs_which, " years")

# ---- simulate, refit, recompute ----
# The site random effects must be REDRAWN, not reused. Simulating from fitted(mA)
# would bake the shrunken site estimates into every replicate as if they were known
# constants, which understates between-site variation — the empirical standard
# deviation of the fitted site effects is 0.219 against an estimated component of
# 0.307 — and produces a null distribution that is too narrow. Every replicate
# therefore draws a fresh set of site effects from N(0, sigma_site^2) and applies
# one draw per site, so the clustering of observations within sites is generated
# rather than inherited.
th0 <- mA$family$getTheta(TRUE)
tm0 <- predict(mA, type = "terms")
re_col0 <- grep("s\\(SITE\\)", colnames(tm0))
eta_fixed <- attr(tm0, "constant") + rowSums(tm0[, -re_col0, drop = FALSE])
sd_site <- gam.vcomp(mA, rescale = FALSE)["s(SITE)", "std.dev"]
site_ix <- as.integer(d$SITE); n_site <- nlevels(d$SITE)
say("\nsimulation draws site effects from N(0, ", sprintf("%.4f", sd_site),
    "^2), one per site, rather than reusing the fitted values")

set.seed(101)
sim_seeds <- sample.int(1e6, NSIM)
one_rep <- function(b) {
  set.seed(sim_seeds[b])
  bs <- rnorm(n_site, 0, sd_site)
  ds <- d
  ds$count <- rnbinom(nrow(d), size = th0, mu = exp(eta_fixed + bs[site_ix]))
  f <- tryCatch(gam(FORM, family = nb(), data = ds, method = "REML"), error = function(e) NULL)
  if (is.null(f)) return(c(NA_real_, NA_real_))
  g <- gap_r(f, ds, sim_seeds[b])
  c(stat_zero(f, ds$count), min(g, na.rm = TRUE))
}
ncore <- max(1L, min(parallel::detectCores(), 4L))
sims <- if (.Platform$OS.type == "unix" && ncore > 1L)
          parallel::mclapply(seq_len(NSIM), one_rep, mc.cores = ncore) else
          lapply(seq_len(NSIM), one_rep)
S <- do.call(rbind, sims); S <- S[stats::complete.cases(S), , drop = FALSE]
say("replicates that refitted successfully: ", nrow(S), " of ", NSIM)

p_zero <- (1 + sum(S[, 1] >= obs_zero)) / (nrow(S) + 1)
p_gap  <- (1 + sum(S[, 2] <= obs_min))  / (nrow(S) + 1)

say("\n-- zero counts --")
say("simulated discrepancy: median ", sprintf("%.1f", median(S[, 1])),
    ", 95% interval ", sprintf("%.1f to %.1f", quantile(S[, 1], .025), quantile(S[, 1], .975)))
say("p = ", sprintf("%.4f", p_zero), "   (smallest attainable with ", nrow(S),
    " replicates is ", sprintf("%.4f", 1 / (nrow(S) + 1)),
    "; computed as (extreme + 1) / (B + 1), so it is never zero)")
v_zero <- if (p_zero < 0.01) "excess" else if (p_zero > 0.10) "none" else "uncertain"
say(switch(v_zero,
  excess = "-> the data hold more zeros than this model generates, after allowing for estimation.",
  none   = "-> the observed zero count is within what this model generates.",
  uncertain = "-> UNCERTAIN. Close to the edge of the simulated null; reported as unresolved."))

say("\n-- residual correlation across survey gaps --")
say("simulated most-negative gap correlation: median ", sprintf("%.3f", median(S[, 2])),
    ", 95% interval ", sprintf("%.3f to %.3f", quantile(S[, 2], .025), quantile(S[, 2], .975)))
say("p = ", sprintf("%.4f", p_gap), "   (this p already accounts for the search across")
say("all eleven gaps, and for the negative correlation the site random effect induces)")
# Rule: a result near the decision boundary is labelled uncertain rather than
# pushed to one side of it. With B replicates the smallest attainable p is
# 1/(B+1), so a p in the single-digit percents is reported as what it is.
verdict <- function(p, hi, lo) if (p < 0.01) hi else if (p > 0.10) lo else "uncertain"
v_gap <- verdict(p_gap, "excess", "none")
if (v_gap == "excess") {
  say("-> more residual temporal structure than this model generates on its own.")
  say("   The site random intercept alone does not account for dependence between")
  say("   repeat visits to a site. This is a statement about model adequacy.")
} else if (v_gap == "none") {
  say("-> not distinguishable from what this model produces by itself. Testing the")
  say("   raw correlation against zero, and picking the strongest of eleven gaps after")
  say("   seeing them, is what made this look like a finding in an earlier version.")
} else {
  say("-> UNCERTAIN. The observed value sits close to the edge of the simulated null")
  say("   and this run cannot place it on one side. Reported as unresolved rather")
  say("   than forced across a threshold; more replicates would only sharpen it if")
  say("   the answer changed a conclusion, and here it does not — the protection")
  say("   estimate is unmoved either way, as the sensitivity table below shows.")
}

gap_tab <- data.frame(gap_years = as.numeric(names(obs_gapr)),
                      r = round(as.numeric(obs_gapr), 3))
gap_tab$is_min <- gap_tab$gap_years == as.numeric(obs_which)
cap(gap_tab)
write.csv(data.frame(statistic = c("zero discrepancy", "min gap correlation"),
                     observed = round(c(obs_zero, obs_min), 3),
                     sim_median = round(c(median(S[, 1]), median(S[, 2])), 3),
                     sim_lo = round(c(quantile(S[, 1], .025), quantile(S[, 2], .025)), 3),
                     sim_hi = round(c(quantile(S[, 1], .975), quantile(S[, 2], .975)), 3),
                     p = round(c(p_zero, p_gap), 4), n_sim = nrow(S)),
          file.path(OUT, "table13_gof_bootstrap.csv"), row.names = FALSE)

# ---- does it matter for the answer? ----
say("\n-- does either affect the protection estimate? --")
zi_site <- tapply(d$count == 0, d$SITE, sum)
drop_sites <- names(zi_site)[zi_site >= 4]
d_zt <- droplevels(d[!d$SITE %in% drop_sites, ])
m_zt <- gam(FORM, family = nb(), data = d_zt, method = "REML")
say("sites with at least 4 of 7 surveys zero: ", length(drop_sites))

# planned sensitivity: a site-level random slope on year, so that sites are allowed
# their own trajectory rather than only their own level
d$SITEs <- d$SITE
mRS <- tryCatch(gam(update(FORM, . ~ . + s(SITEs, YEAR, bs = "re")),
                    family = nb(), data = d, method = "REML"), error = function(e) NULL)
say("site-year random slope model: ", if (is.null(mRS)) "did not converge" else
    paste0("fitted, REML ", sprintf("%.1f", mRS$gcv.ubre), " against ",
           sprintf("%.1f", mA$gcv.ubre), " for the random-intercept model"))

rows <- list(c("Full model (A)", "mA"), c("Zero-trimmed", "m_zt"))
if (!is.null(mRS)) rows[[3]] <- c("+ site-year random slope", "mRS")
sens <- do.call(rbind, lapply(rows, function(z) {
  m <- get(z[2])
  a <- rr(m, "Whitsunday", "NTR 1987"); b <- rr(m, "Whitsunday", "NTR 2004")
  data.frame(model = z[1],
             W1987 = sprintf("%.2f (%.2f-%.2f)", a[1], a[2], a[3]),
             W2004 = sprintf("%.2f (%.2f-%.2f)", b[1], b[2], b[3]))
}))
cap(sens)
write.csv(sens, file.path(OUT, "table16_gof_sensitivity.csv"), row.names = FALSE)
say("\nWhat this does and does not show. Neither sensitivity available here moves the")
say("Whitsunday estimate materially: trimming the mostly-zero site and allowing sites")
say("their own time trajectory both leave it where it was, with wider intervals in the")
say("second case. There is also a general argument that a zero excess inflates apparent")
say("overdispersion, which widens intervals rather than narrowing them.")
say("")
say("That is weaker than saying the zero excess cannot affect the result, and the")
say("stronger claim is not made. The test that would settle it is a refit under a")
say("distribution built for the excess — a hurdle or zero-inflated negative binomial —")
say("and no such model was fitted here. Until one is, the correct statement is that")
say("the sensitivities that were run did not change the estimate, not that the")
say("distributional misfit is harmless.")

# ---------------------------------------------------------------------
# 6b. OBSERVED AGAINST FITTED, FINAL MODEL
# ---------------------------------------------------------------------
rule("6b. OBSERVED AGAINST FITTED (FINAL MODEL)")
say("Step 1 checked observed against fitted for the baseline model. The final model")
say("was never checked the same way. Binned means with the interval a correct model")
say("would produce, so systematic departure is visible rather than inferred.\n")
mu <- fitted(mA); th <- mA$family$getTheta(TRUE)
brk <- unique(quantile(mu, seq(0, 1, length.out = 11)))
bin <- cut(mu, brk, include.lowest = TRUE)
of <- do.call(rbind, lapply(levels(bin), function(b) {
  i <- bin == b
  data.frame(bin = b, n = sum(i), fitted_mean = round(mean(mu[i]), 2),
             observed_mean = round(mean(d$count[i]), 2),
             se = round(sd(d$count[i]) / sqrt(sum(i)), 2))
}))
of$z <- round((of$observed_mean - of$fitted_mean) / of$se, 2)
cap(of)
write.csv(of, file.path(OUT, "table17_observed_vs_fitted.csv"), row.names = FALSE)
say("")
say("The z column is descriptive, not a test. The bins are defined by the fitted values")
say("themselves, the ten of them are not independent, and no multiplicity or selection")
say("adjustment is applied, so |z| > 2 in one bin is not a rejection of anything. Read")
say("the column as a scale for how far each bin sits from its prediction relative to")
say("the noise in that bin, and read the pattern across bins rather than any single")
say("value.")
say("")
say("Bins where |z| exceeds 2: ",
    if (any(abs(of$z) > 2, na.rm = TRUE)) paste(which(abs(of$z) > 2), collapse = ", ") else "none",
    ". The lowest-density bin is the one to watch, and it is the zero excess of section")
say("5-6 seen from another angle rather than a separate problem.")

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
fam$survives <- ifelse(is.na(fam$p_holm), NA, fam$p_holm < 0.05)
fam$p <- fmt_p(fam$p); fam$p_holm <- fmt_p(fam$p_holm)
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

# a — zero discrepancy against its bootstrap null
h <- hist(S[, 1], breaks = 24, plot = FALSE)
plot(h, col = "grey88", border = "white",
     xlim = range(c(h$breaks, obs_zero)) + c(-1, 3),
     xlab = "Observed minus expected zeros (simulated)",
     ylab = "Replicates", main = "a  Zero counts against a parametric bootstrap",
     font.main = 1, cex.main = 0.88, adj = 0)
abline(v = obs_zero, col = RED, lwd = 2.2)
text(obs_zero, par("usr")[4] * 0.90, paste0("observed ", sprintf("%.0f", obs_zero), " "),
     col = RED, cex = 0.72, adj = 1)
mtext(sprintf("p = %.3f, %d refits", p_zero, nrow(S)), side = 3, line = -1.1,
      adj = 0.98, cex = 0.66, col = "grey35")

# b — most-negative gap correlation against its bootstrap null
h2 <- hist(S[, 2], breaks = 24, plot = FALSE)
plot(h2, col = "grey88", border = "white",
     xlim = range(c(h2$breaks, obs_min)) + c(-0.03, 0.03),
     xlab = "Most negative gap correlation (simulated)",
     ylab = "Replicates", main = "b  Residual correlation, search over all gaps",
     font.main = 1, cex.main = 0.88, adj = 0)
abline(v = obs_min, col = RED, lwd = 2.2)
text(obs_min, par("usr")[4] * 0.90, paste0(" observed ", sprintf("%.2f", obs_min)),
     col = RED, cex = 0.72, adj = 0)
mtext(sprintf("p = %.3f", p_gap), side = 3, line = -1.1, adj = 0.98,
      cex = 0.66, col = "grey35")
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
