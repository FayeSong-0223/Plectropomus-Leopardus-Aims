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
rule("4. HOW MUCH DOES EACH BLOCK CONTRIBUTE?")
say("Q2 asks how much of the variation the covariates explain. An earlier version of")
say("this analysis answered it by decomposing the variance of the linear predictor")
say("into per-block shares. That answer has been withdrawn, for three reasons found")
say("on audit and recorded in DECISIONS.md:")
say("")
say("  1. It was not invariant to the coding of the factor contrasts. Refitting the")
say("     identical model under sum-to-zero rather than treatment contrasts moved the")
say("     management share from 0.481 to 0.391 and space/time from 0.121 to 0.345,")
say("     with identical fitted values and identical log-likelihood. A quantity that")
say("     moves when nothing about the fit moves is a property of the parameterisation,")
say("     not of the data.")
say("  2. It partitioned only the systematic part. The linear predictor has variance")
say("     0.533; the negative binomial observation process contributes roughly 0.455")
say("     more on the same scale. Around half of the variation sat outside the")
say("     partition and the shares were silently conditional on that.")
say("  3. The site block used the empirical variance of the shrunken random-effect")
say("     predictions, 0.048, against an estimated variance component of 0.094. Site")
say("     was understated by about half, and every other block inflated against it.")
say("")
say("What replaces it is a drop-one-block comparison: refit the model without each")
say("block in turn and record how far deviance explained falls, in percentage points.")
say("This is invariant to contrast coding, it is a statement about the fitted model")
say("rather than about total ecological variation, and it can carry an interval.")
say("It is not additive: blocks share variation, so the drops do not sum to the")
say("model's deviance explained and should never be presented as shares of it.\n")

blocks <- list(
  "Management (H1/Q3)"        = ". ~ . - REGION:NTR - NTR",
  "Habitat (H2)"              = ". ~ . - s(rugosity, k = 5) - s(LHC, k = 5) - depth",
  "Environment (H3)"          = ". ~ . - s(kd490, k = 5) - s(maxDHW, k = 5) - s(Cyclone, k = 5)",
  "Regional year trends"      = ". ~ . - s(YEAR, by = REGION, k = 5)",
  "Site random effect (H4)"   = ". ~ . - s(SITE, bs = \"re\")",
  "Wave exposure (confounder)"= ". ~ . - EXPOSURE")

fA   <- formula(mA)
fNoRE <- update(fA, . ~ . - s(SITE, bs = "re"))

gsafe <- function(fm, dat) tryCatch(gam(fm, family = nb(), data = dat, method = "REML"),
                                    error = function(e) NULL)
LABS <- c(names(blocks), "Management, no site random effect")
dev_drops <- function(dat) {
  out <- rep(NA_real_, length(LABS))
  full <- gsafe(fA, dat)
  if (!is.null(full)) {
    de0 <- summary(full)$dev.expl
    for (i in seq_along(blocks)) {
      m <- gsafe(update(fA, as.formula(blocks[[i]])), dat)
      if (!is.null(m)) out[i] <- (de0 - summary(m)$dev.expl) * 100
    }
  }
  # Management is constant within a site, so with a site random effect in the
  # model the two compete for the same between-site variation and dropping
  # management simply hands its work to the random effect. To measure what
  # management explains at all, it has to be dropped from a model where nothing
  # else is absorbing site identity.
  f0 <- gsafe(fNoRE, dat)
  f1 <- gsafe(update(fNoRE, . ~ . - REGION:NTR - NTR), dat)
  if (!is.null(f0) && !is.null(f1))
    out[length(LABS)] <- (summary(f0)$dev.expl - summary(f1)$dev.expl) * 100
  setNames(out, LABS)
}

obs <- dev_drops(d)
say("deviance explained, full model: ", sprintf("%.1f%%", summary(mA)$dev.expl * 100))
say("deviance explained, same model without the site random effect: ",
    sprintf("%.1f%%", summary(gam(fNoRE, family = nb(), data = d, method = "REML"))$dev.expl * 100))

# Cluster bootstrap: sites are resampled with replacement, because observations
# within a site are not independent. Resampled copies of the same site get
# distinct labels so the random effect treats them as separate draws. The
# resample indices are drawn up front under a fixed seed, so the result does not
# depend on how the work is later divided between cores.
NBOOT <- 200
sites <- levels(d$SITE)
set.seed(2)
boot_sets <- lapply(seq_len(NBOOT), function(b) sample(sites, length(sites), replace = TRUE))

make_boot <- function(sel) {
  parts <- lapply(seq_along(sel), function(j) {
    z <- d[d$SITE == sel[j], ]; z$SITE <- paste0(sel[j], "_", j); z
  })
  z <- do.call(rbind, parts); z$SITE <- factor(z$SITE); z
}

say("\nbootstrapping ", NBOOT, " site-level resamples (", length(blocks) + 3,
    " model fits each) — this is by far the slow part of the pipeline,")
say("of the order of half an hour on two cores and a few minutes on eight")
run_one <- function(b) dev_drops(make_boot(boot_sets[[b]]))
ncore <- max(1L, min(parallel::detectCores(), 4L))
bt <- if (.Platform$OS.type == "unix" && ncore > 1L)
        parallel::mclapply(seq_len(NBOOT), run_one, mc.cores = ncore) else
        lapply(seq_len(NBOOT), run_one)
bt <- do.call(rbind, lapply(bt, function(z)
  if (is.numeric(z) && length(z) == length(LABS)) z else rep(NA_real_, length(LABS))))
ok <- rowSums(is.na(bt)) == 0
say("replicates that converged for every model: ", sum(ok), " of ", NBOOT)

ci <- t(apply(bt[ok, , drop = FALSE], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
dd <- data.frame(block = LABS,
                 dev_expl_drop_pp = round(as.numeric(obs), 2),
                 lo = round(ci[, 1], 2), hi = round(ci[, 2], 2))
dd$excludes_zero <- !(dd$lo < 0 & dd$hi > 0)
dd <- dd[order(-dd$dev_expl_drop_pp), ]
say("\n-- fall in deviance explained when each block is removed (percentage points) --")
cap(dd)
write.csv(dd, file.path(OUT, "table5_block_contributions.csv"), row.names = FALSE)

say("\nThese are conditional, non-additive measures of each block's contribution to")
say("this model. They are not shares of total ecological variation and do not sum")
say("to anything meaningful. A block whose interval includes zero is not")
say("distinguishable from a block that contributes nothing.")
say("")
say("The management row is the one that repays reading carefully, and it is the")
say("honest answer to Q2. Protection is fixed for the whole history of a site, so")
say("in a model that already gives every site its own intercept, management and the")
say("random effect are competing to explain the same between-site differences.")
say("Remove management and the random effect simply takes the work over, which is")
say("why its drop sits at or below zero. That is not evidence of no effect — Step 3")
say("shows the Whitsunday estimate is stable across every specification — it means")
say("deviance explained cannot separate the two, and any method that appears to")
say("separate them is reporting an artefact of how the model was parameterised.")
say("The last row measures management where it is identifiable: in a model with no")
say("site random effect at all.")
say("")
say("H4 said unexplained site-level variation would be large relative to the")
say("measured covariates. On these numbers it is: the site random effect is the")
say("largest single contributor, and no covariate block approaches it.")

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
  plot(mm$m, select = k, shade = TRUE, shade.col = adjustcolor(TEAL, 0.20),
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
