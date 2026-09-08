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
say("Effective degrees of freedom are reported so over-fitting is visible.")
say("")
say("That justification was later tested rather than left as an assertion: Step 5")
say("refits at k = 10 and k = 20 for kd490 and at k = 6 for the year smooths, the")
say("largest the survey years allow. The basis warnings do not clear at any of them")
say("and the protection estimate does not move, so the residual pattern they detect")
say("is not something a larger basis represents. See table 18.\n")

mA <- gam(count ~ REGION * NTR + EXPOSURE + depth +
            s(YEAR, by = REGION, k = 5) +
            s(rugosity, k = 5) + s(LHC, k = 5) +
            s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) +
            s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")

mB <- gam(count ~ REGION * NTR + EXPOSURE + depth + RY +
            s(rugosity, k = 5) + s(LHC, k = 5) +
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
say("Five models are compared, from the Step 1 baseline to the fully adjusted fits.")
say("")
say("The first rung is the Step 1 model exactly as Step 1 fitted it, with no depth")
say("term, so that its numbers reproduce table 3 rather than merely resembling them.")
say("Depth enters as its own rung, because it moves the Whitsunday estimate more")
say("than any other single covariate and that should be visible rather than folded")
say("into a baseline.\n")

mS1 <- gam(count ~ REGION * NTR + EXPOSURE + factor(YEAR) + s(SITE, bs = "re"),
           family = nb(), data = d, method = "REML")
m0 <- gam(count ~ REGION * NTR + EXPOSURE + depth + factor(YEAR) + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")
mH <- gam(count ~ REGION * NTR + EXPOSURE + depth + factor(YEAR) +
            s(rugosity, k = 5) + s(LHC, k = 5) + s(SITE, bs = "re"),
          family = nb(), data = d, method = "REML")

# Intervals use the unconditional covariance matrix throughout: it adds the
# uncertainty in the smoothing parameters, which the conditional matrix treats
# as known. For protection the difference is negligible (these are parametric
# contrasts); for anything read off a smooth it is not, and using one rule
# everywhere avoids two standards in one report.
rr <- function(m, region, level) {
  b <- coef(m); V <- vcov(m, unconditional = TRUE); nm <- names(b); k <- rep(0, length(b))
  k[match(paste0("NTR", level), nm)] <- 1
  ix <- paste0("REGIONWhitsunday:NTR", level)
  if (region == "Whitsunday" && ix %in% nm) k[match(ix, nm)] <- 1
  e <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(exp(e), exp(e - 1.96 * se), exp(e + 1.96 * se))
}
mods <- list("1 Step 1 baseline (no depth)" = mS1, "2 + depth" = m0, "3 + habitat" = mH,
             "4 + environment, smooth year (A)" = mA, "5 region-year saturated (B)" = mB)
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
rule("4. HOW MUCH DOES EACH BLOCK CONTRIBUTE?  (Q2)")
say("Q2 asks how much of the variation the covariates explain. Two earlier answers")
say("have been withdrawn, and the reasons matter more than the numbers did.")
say("")
say("Withdrawn 1: a link-scale variance decomposition reporting management at 0.481.")
say("  It was not invariant to the coding of the factor contrasts. Refitting the")
say("  identical model under sum-to-zero moved management to 0.391 and space/time")
say("  from 0.121 to 0.345, with identical fitted values and log-likelihood. It also")
say("  partitioned only the systematic half of the variation, and used the shrunken")
say("  random-effect predictions for the site block, understating it by about half.")
say("")
say("Withdrawn 2: a drop-one-block comparison in which management was measured by")
say("  removing it from a model with no site random effect. That number, 14.09")
say("  percentage points, is not a management contribution. With no random effect")
say("  present, the management terms absorb every between-site difference that")
say("  happens to align with zoning — habitat, history, location, anything fixed")
say("  about a site. It measures confounding as much as management, and reporting")
say("  it as management's share would have repeated the error it was meant to fix.")
say("")
say("What is reported instead is set out below. It begins with the thing that has to")
say("be said plainly rather than measured.\n")

# ---------------------------------------------------------------------
say(strrep("-", 70))
say("4a. WHAT IS NOT IDENTIFIABLE, AND WHY")
say(strrep("-", 70))
nvary <- sum(tapply(as.integer(d$NTR), d$SITE, function(z) length(unique(z))) > 1)
say("Protection changes within ", nvary, " of ", nlevels(d$SITE), " sites.")
say("")
say("That single fact settles it. Protection is a fixed property of a site across the")
say("whole series, so its entire contribution is between-site. A site random intercept")
say("is also, by construction, entirely between-site. The two are estimated from the")
say("same 71 degrees of freedom, and no partition of variance or of deviance can say")
say("how much of a between-site difference belongs to zoning rather than to whatever")
say("else distinguishes those sites. Drop management and the random effect absorbs it;")
say("drop the random effect and management absorbs everything else.")
say("")
say("This is a limit of the design, not of the method, and it cannot be fixed by")
say("choosing a better statistic. What CAN be estimated is the size of the protection")
say("contrast, with an interval, which is Step 1 and Step 3's business — and how much")
say("the model's ability to predict a site it has never seen depends on knowing that")
say("site's zoning. That second question is answerable, and it is asked in 4c.")

# ---------------------------------------------------------------------
say("\n", strrep("-", 70))
say("4b. DROP-ONE-BLOCK AT FIXED DISPERSION")
say(strrep("-", 70))
say("Each block is removed in turn and the fall in deviance explained recorded.")
say("")
say("One correction to how this was done before: nb() re-estimates theta for every")
say("model, so a reduced model can absorb its lost structure into a smaller theta and")
say("the two deviances are then measured on different scales, which makes the")
say("difference uninterpretable. Every model below is fitted with theta held at the")
say("full model's estimate, using negbin(), so all deviances share one scale.")
say("")
say("Management is included for completeness and is expected to be near zero for the")
say("reason given in 4a. It is not evidence of no effect.\n")

TH <- mA$family$getTheta(TRUE)
say("theta from the full model, held fixed throughout: ", sprintf("%.4f", TH))

blocks <- list(
  "Management (H1/Q3)"        = ". ~ . - REGION:NTR - NTR",
  "Habitat (H2)"              = ". ~ . - s(rugosity, k = 5) - s(LHC, k = 5) - depth",
  "Environment (H3)"          = ". ~ . - s(kd490, k = 5) - s(maxDHW, k = 5) - s(Cyclone, k = 5)",
  "Regional year trends"      = ". ~ . - s(YEAR, by = REGION, k = 5)",
  "Site random effect (H4)"   = ". ~ . - s(SITE, bs = \"re\")",
  "Wave exposure (confounder)"= ". ~ . - EXPOSURE")

fA <- formula(mA)
gfix <- function(fm, dat) tryCatch(gam(fm, family = negbin(TH), data = dat, method = "REML"),
                                   error = function(e) NULL)
dev_drops <- function(dat) {
  out <- rep(NA_real_, length(blocks))
  full <- gfix(fA, dat)
  if (!is.null(full)) {
    de0 <- summary(full)$dev.expl
    for (i in seq_along(blocks)) {
      m <- gfix(update(fA, as.formula(blocks[[i]])), dat)
      if (!is.null(m)) out[i] <- (de0 - summary(m)$dev.expl) * 100
    }
  }
  setNames(out, names(blocks))
}
obs <- dev_drops(d)
say("deviance explained at fixed theta, full model: ",
    sprintf("%.1f%%", summary(gfix(fA, d))$dev.expl * 100))

NBOOT <- 200
sites <- levels(d$SITE)
set.seed(2)
boot_sets <- lapply(seq_len(NBOOT), function(b) sample(sites, length(sites), replace = TRUE))
make_boot <- function(sel) {
  z <- do.call(rbind, lapply(seq_along(sel), function(j) {
    w <- d[d$SITE == sel[j], ]; w$SITE <- paste0(sel[j], "_", j); w }))
  z$SITE <- factor(z$SITE); z
}
say("\nbootstrapping ", NBOOT, " site-level resamples at ", length(blocks) + 1,
    " fits each; resample indices are drawn up front under a fixed seed, so the")
say("result does not depend on how the work is divided between cores")
ncore <- max(1L, min(parallel::detectCores(), 4L))
bt <- if (.Platform$OS.type == "unix" && ncore > 1L)
        parallel::mclapply(seq_len(NBOOT), function(b) dev_drops(make_boot(boot_sets[[b]])),
                           mc.cores = ncore) else
        lapply(seq_len(NBOOT), function(b) dev_drops(make_boot(boot_sets[[b]])))
bt <- do.call(rbind, lapply(bt, function(z)
  if (is.numeric(z) && length(z) == length(blocks)) z else rep(NA_real_, length(blocks))))
ok <- rowSums(is.na(bt)) == 0
say("replicates that converged for every model: ", sum(ok), " of ", NBOOT)

ci <- t(apply(bt[ok, , drop = FALSE], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
dd <- data.frame(block = names(blocks), dev_expl_drop_pp = round(as.numeric(obs), 2),
                 lo = round(ci[, 1], 2), hi = round(ci[, 2], 2))
dd$excludes_zero <- !(dd$lo < 0 & dd$hi > 0)
dd <- dd[order(-dd$dev_expl_drop_pp), ]; rownames(dd) <- NULL
say("\n-- fall in deviance explained when each block is removed (percentage points) --")
cap(dd)
write.csv(dd, file.path(OUT, "table5_block_contributions.csv"), row.names = FALSE)
say("\nConditional and non-additive. Blocks share variation, the drops do not sum to")
say("anything, and they are not shares of total variation.")

# ---------------------------------------------------------------------
say("\n", strrep("-", 70))
say("4c. OUT-OF-SAMPLE: WHAT DOES KNOWING A NEW SITE'S ZONING BUY?")
say(strrep("-", 70))
say("The in-sample comparisons above cannot separate management from site identity")
say("because both are fitted to the same sites. Prediction to a site the model has")
say("never seen removes the need to ESTIMATE that site's intercept: it is unknown")
say("and integrated out rather than fitted alongside management. That is a different")
say("thing from removing confounding — whatever distinguished the sites selected for")
say("zoning is still present in the held-out sites, and the comparison below cannot")
say("see it. What follows is about prediction, not about identification.")
say("")
say("Ten folds, split by SITE so that a site is never in both training and test.")
say("Each model is fitted on the training sites and scored on the held-out sites by")
say("the MARGINAL negative binomial log predictive density of individual site-year")
say("counts — the site effect integrated over N(0, sigma^2) rather than set to zero,")
say("since under a log link the density at the median of that distribution is not its")
say("mean. Each model uses its own theta: a predictive comparison is between complete")
say("distributions, so re-estimating theta is correct here.")
say("")
say("This is exploratory. It measures predictive contribution, not a share of")
say("variance, and it does not license a causal reading: a zoning label may predict a")
say("new site well because of what protection does, or because of what protection is")
say("correlated with.\n")

K <- 10
set.seed(7)
fold_of_site <- setNames(sample(rep_len(1:K, length(sites))), sites)
d$.fold <- fold_of_site[as.character(d$SITE)]

# ---- Gauss-Hermite nodes and weights, base R only (Golub-Welsch) -----------
# For an unseen site the random intercept is unknown, so the predictive density is
#   p(y) = INT NB(y; exp(eta + b), theta) * N(b; 0, sigma^2) db
# rather than the density at b = 0. Under a log link those differ: the density at
# the median of the random-effect distribution is not its mean. An earlier version
# evaluated at b = 0 and is corrected here.
#
# Substituting b = sqrt(2) * sigma * x turns the integral into the Gauss-Hermite
# form INT f(x) exp(-x^2) dx, approximated by sum_i w_i f(x_i), so
#   log p(y) = -0.5*log(pi) + logsumexp_i [ log(w_i) + log NB(y; exp(eta + sqrt(2) sigma x_i), theta) ]
# The logsumexp keeps this stable: the individual densities underflow to zero in
# the tails, and summing them on the natural scale would silently lose them.
gauss_hermite <- function(n) {
  i <- seq_len(n - 1); J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i / 2); J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), w = rev(sqrt(pi) * e$vectors[1, ]^2))
}
GH <- gauss_hermite(20)
stopifnot(abs(sum(GH$w) - sqrt(pi)) < 1e-10)          # weights integrate exp(-x^2)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

# The site random effect is deliberately NOT among the blocks compared here. It is
# already integrated out of every prediction, so a "drop the random effect" row would
# not measure site identity's predictive contribution and would invite exactly the
# comparison this section cannot support — management against site identity. That
# comparison is not available from this design by any route (4a).
cv_blocks <- blocks[names(blocks) != "Site random effect (H4)"]
cv_models <- c(list("Full model" = ". ~ ."), cv_blocks)
lpd <- matrix(NA_real_, nrow = nrow(d), ncol = length(cv_models),
              dimnames = list(NULL, names(cv_models)))
sig_fold <- rep(NA_real_, K)

for (k in 1:K) {
  tr <- droplevels(d[d$.fold != k, ]); te <- d[d$.fold == k, ]
  for (mi in seq_along(cv_models)) {
    fm <- if (names(cv_models)[mi] == "Full model") fA else update(fA, as.formula(cv_models[[mi]]))
    fit <- tryCatch(gam(fm, family = nb(), data = tr, method = "REML"), error = function(e) NULL)
    if (is.null(fit)) next
    trm <- predict(fit, newdata = transform(te, SITE = tr$SITE[1]), type = "terms")
    re_col <- grep("s\\(SITE\\)", colnames(trm))
    eta <- attr(trm, "constant") +
           rowSums(trm[, setdiff(seq_len(ncol(trm)), re_col), drop = FALSE])
    th_k <- fit$family$getTheta(TRUE)
    # sigma is re-estimated inside each training fold, from that fold's own fit,
    # so nothing from the held-out sites leaks into the prediction.
    vc <- gam.vcomp(fit, rescale = FALSE)
    sg <- if ("s(SITE)" %in% rownames(vc)) vc["s(SITE)", "std.dev"] else 0
    if (names(cv_models)[mi] == "Full model") sig_fold[k] <- sg
    lpd[d$.fold == k, mi] <- vapply(seq_along(eta), function(r)
      -0.5 * log(pi) + lse(log(GH$w) +
        dnbinom(te$count[r], size = th_k,
                mu = exp(eta[r] + sqrt(2) * sg * GH$x), log = TRUE)),
      numeric(1))
  }
}
say("site-effect SD estimated separately in each training fold: ",
    paste(sprintf("%.3f", sig_fold), collapse = ", "))

base_lpd <- lpd[, "Full model"]
site_of <- as.character(d$SITE)
cvres <- do.call(rbind, lapply(names(cv_blocks), function(nm) {
  dif <- base_lpd - lpd[, nm]
  keep <- is.finite(dif)
  # Clustered standard error. Observations are not independent within a site, so the
  # difference is summed within each site first and the spread taken across the 71
  # site totals. The naive version, which treats all 467 observations as independent,
  # is reported alongside because an earlier version quoted it as though it were the
  # standard error, and the gap between them is the point.
  per_site <- tapply(dif[keep], site_of[keep], sum)
  ns <- length(per_site)
  data.frame(block = nm, delta_lpd = round(sum(dif[keep]), 1),
             se_site = round(sqrt(ns) * sd(per_site), 1),
             se_naive = round(sqrt(sum(keep)) * sd(dif[keep]), 1),
             n_sites = ns, per_obs = round(mean(dif[keep]), 4))
}))
cvres$z_site  <- round(cvres$delta_lpd / cvres$se_site, 2)
cvres$z_naive <- round(cvres$delta_lpd / cvres$se_naive, 2)
cvres <- cvres[order(-cvres$delta_lpd), ]; rownames(cvres) <- NULL
say("\nTotal log predictive density of the full model across all held-out sites: ",
    sprintf("%.1f", sum(base_lpd[is.finite(base_lpd)])))
say("\n-- loss in log predictive density when each block is removed --")
say("   positive = removing that block makes prediction to an unseen site worse.")
say("   se_site clusters on the ", cvres$n_sites[1], " sites; se_naive treats all")
say("   observations as independent and is shown only for comparison.\n")
cap(cvres)
write.csv(cvres, file.path(OUT, "table15_cv_block_contributions.csv"), row.names = FALSE)

say("\nHow to read this, and how not to.")
say("  It is a predictive statement: knowing a new site's zoning improves the")
say("  prediction of its counts by this much on the log scale. It does not partition")
say("  variance, the entries do not sum to anything, and a large value is not evidence")
say("  that management contributes more variance than site identity — that comparison")
say("  is unavailable in this design by any route, as 4a sets out.")
say("  It is not causal. Zoning may predict an unseen site well because of what")
say("  protection does, or because of whatever the zoning process selected for. This")
say("  analysis cannot tell those apart, and nothing here should be read as trying to.")

# ---------------------------------------------------------------------
# 5. FIGURES
# ---------------------------------------------------------------------
rule("5. FIGURES")

## Figure 6 — partial effects
png(file.path(OUT, "fig6_partial_effects.png"), width = 2100, height = 1400, res = 220)
par(mfrow = c(2, 3), mar = c(3.6, 3.8, 2.4, 0.8)); base_par()
sm <- c("s(rugosity)", "s(LHC)", "s(kd490)", "s(maxDHW)", "s(Cyclone)")
lab <- c("Rugosity index", "Live hard coral (%)",
         "kd490 (turbidity)", "Max degree heating weeks", "Cyclone exposure index")
st <- summary(mA)$s.table
for (i in seq_along(sm)) {
  j <- match(sm[i], rownames(st))
  plot(mA, select = which(sapply(mA$smooth, function(z) z$label) == sm[i]),
       unconditional = TRUE,   # bands must match the intervals reported in Step 4
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

## Figure 7 — block contributions, with bootstrap intervals
png(file.path(OUT, "fig7_block_contributions.png"), width = 1900, height = 1150, res = 220)
par(mar = c(4.2, 11.0, 2.8, 1.4)); base_par()
o  <- dd[order(dd$dev_expl_drop_pp), ]
xr <- range(c(0, o$lo, o$hi)) * c(1, 1.10)
bp <- barplot(o$dev_expl_drop_pp, horiz = TRUE, names.arg = o$block, col = TEAL,
              border = NA, xlim = xr, las = 1, cex.names = 0.78,
              xlab = "Fall in deviance explained when removed (percentage points)",
              main = "Contribution of each block to this model", font.main = 1,
              cex.main = 0.98, adj = 0)
arrows(o$lo, bp, o$hi, bp, angle = 90, code = 3, length = 0.03, col = "grey30", xpd = NA)
abline(v = 0, col = "grey60", lty = 2)
mtext("Conditional and non-additive: these are not shares of total variation",
      side = 3, line = -0.1, adj = 0, cex = 0.72, col = "grey35")
dev.off(); say("wrote fig7_block_contributions.png")

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
  axis(1, at = xs, labels = c("Step 1", "+depth", "+hab", "+env", "sat"), cex.axis = 0.75)
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
mA2 <- gam(count ~ REGION * NTR + EXPOSURE + depth +
             s(YEAR, by = REGION, k = 5) +
             s(rugosity, k = 5) + s(LHC, k = 5) +
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
  plot(mm$m, select = k, unconditional = TRUE,
       shade = TRUE, shade.col = adjustcolor(TEAL, 0.20),
       col = TEAL, lwd = 2, rug = TRUE, ylim = c(-1.2, 1.2),
       xlab = "Max degree heating weeks", ylab = "Effect on log density",
       main = mm$t, font.main = 1, cex.main = 1.0, adj = 0)
  abline(h = 0, col = "grey60", lty = 2)
}
dev.off(); say("\nwrote fig9_thermal_sensitivity.png")

say("\nRead this comparison before quoting any thermal effect. If the curve changes")
say("shape when one survey is removed, the effect is that survey, not the covariate.")

# Record the environment. mgcv ships with R but its version tracks the R
# version, and REML fitting and the nb() family have both changed across
# releases, so "no packages required" is not the same as "no versions to
# reconcile". Anyone reproducing these numbers needs to know what produced them.
say("\n", strrep("-", 70))
say("environment: ", R.version.string, " | mgcv ", as.character(packageVersion("mgcv")),
    " | platform ", R.version$platform)
say(strrep("-", 70))

writeLines(log_lines, file.path(OUT, "stage2_step2_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
