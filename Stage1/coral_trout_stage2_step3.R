# =====================================================================
#  Plectropomus leopardus — Stage 2, Step 3: sensitivity and influence
#  Palm and Whitsunday inshore reefs, Great Barrier Reef, 2007-2018
#
#  Ziqi (Faye) Song
#
#  Step 2 produced two headline results — a protection effect around 3x
#  in Whitsunday and none in Palm, and a variance decomposition in which
#  management dominates — but also showed the model is sensitive to a
#  single survey (Whitsunday 2017). This step interrogates both.
#
#   1. Influence   every region-year dropped in turn, 13 refits, to see
#                  which surveys the estimates actually depend on.
#   2. Restriction the sheltered stratum alone, where protection and
#                  exposure are far better balanced. Removes the confound
#                  without assuming a functional form for it.
#   3. Separation  each region fitted on its own.
#   4. Cyclone     a lagged term, since disturbance effects need not be
#                  contemporaneous.
#   5. Family      Tweedie on the density scale, to confirm nothing rests
#                  on the negative binomial choice.
#
#  Dependencies: mgcv only.
#  Run:  Rscript coral_trout_stage2_step3.R
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
TEAL <- "#2E7D6F"; RED <- "#C1655A"; BLUE <- "#4A6FA5"; GREY <- "grey55"

# ---------------------------------------------------------------------
# 1. DATA
# ---------------------------------------------------------------------
rule("1. DATA")
f <- file.path(OUT, "analysis_dataset.csv")
if (!file.exists(f)) stop("outputs/analysis_dataset.csv not found. Run coral_trout_stage1.R first.")
d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
names(d)[names(d) == "LHC_%"] <- "LHC"
names(d)[names(d) == "LCC_%"] <- "LCC"
names(d)[names(d) == "Corrected depth"] <- "depth"
d$REGION   <- factor(d$REGION,   levels = c("Palm", "Whitsunday"))
d$NTR      <- factor(d$NTR,      levels = c("Fished", "NTR 1987", "NTR 2004"))
d$EXPOSURE <- factor(d$EXPOSURE, levels = c("Sheltered", "Semi-Exposed", "Exposed"))
d$SITE     <- factor(d$SITE); d$RY <- factor(paste(d$REGION, d$YEAR))
say("rows ", nrow(d), " | sites ", nlevels(d$SITE), " | region-years ", nlevels(d$RY))

# lagged cyclone: the previous survey's value at the same site
d <- d[order(d$SITE, d$YEAR), ]
d$Cyclone_lag <- ave(d$Cyclone, d$SITE, FUN = function(z) c(NA, z[-length(z)]))
say("lagged cyclone constructed (previous survey at the same site); ",
    sum(is.na(d$Cyclone_lag)), " first-visit rows have no lag")

FORM <- count ~ REGION * NTR + EXPOSURE + s(YEAR, by = REGION, k = 5) +
  s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
  s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) + s(SITE, bs = "re")

fitA <- function(dat, form = FORM)
  tryCatch(gam(form, family = nb(), data = dat, method = "REML"), error = function(e) NULL)

rr <- function(m, region, level) {
  if (is.null(m)) return(c(NA, NA, NA))
  b <- coef(m); V <- vcov(m); nm <- names(b); k <- rep(0, length(b))
  main <- paste0("NTR", level); if (!main %in% nm) return(c(NA, NA, NA))
  k[match(main, nm)] <- 1
  ix <- paste0("REGIONWhitsunday:NTR", level)
  if (region == "Whitsunday") { if (!ix %in% nm) return(c(NA, NA, NA)); k[match(ix, nm)] <- 1 }
  e <- sum(k * b); se <- sqrt(as.numeric(t(k) %*% V %*% k))
  c(exp(e), exp(e - 1.96 * se), exp(e + 1.96 * se))
}
sm_stat <- function(m, lab, what = "p-value") {
  if (is.null(m)) return(NA)
  st <- summary(m)$s.table; j <- match(lab, rownames(st)); if (is.na(j)) return(NA)
  st[j, what]
}

mFull <- fitA(d)
say("full-data model refitted for reference; deviance explained ",
    sprintf("%.1f%%", summary(mFull)$dev.expl * 100))

# ---------------------------------------------------------------------
# 2. INFLUENCE: DROP EACH REGION-YEAR IN TURN
# ---------------------------------------------------------------------
rule("2. INFLUENCE — LEAVE ONE REGION-YEAR OUT")
say("Thirteen refits, each omitting one region-year. If a conclusion depends on a")
say("single survey, it appears here as one point far from the rest.\n")

lvls <- levels(d$RY)
inf <- do.call(rbind, lapply(lvls, function(L) {
  m <- fitA(droplevels(d[d$RY != L, ]))
  data.frame(dropped = L, n = sum(d$RY != L),
             W1987 = rr(m, "Whitsunday", "NTR 1987")[1],
             W2004 = rr(m, "Whitsunday", "NTR 2004")[1],
             P1987 = rr(m, "Palm", "NTR 1987")[1],
             dhw_p = sm_stat(m, "s(maxDHW)"),
             kd_p  = sm_stat(m, "s(kd490)"),
             rug_p = sm_stat(m, "s(rugosity)"))
}))
full <- data.frame(dropped = "(none — full data)", n = nrow(d),
                   W1987 = rr(mFull, "Whitsunday", "NTR 1987")[1],
                   W2004 = rr(mFull, "Whitsunday", "NTR 2004")[1],
                   P1987 = rr(mFull, "Palm", "NTR 1987")[1],
                   dhw_p = sm_stat(mFull, "s(maxDHW)"),
                   kd_p  = sm_stat(mFull, "s(kd490)"),
                   rug_p = sm_stat(mFull, "s(rugosity)"))
inf_all <- rbind(full, inf)
inf_all[, 3:8] <- round(inf_all[, 3:8], 3)
cap(inf_all)
write.csv(inf_all, file.path(OUT, "table6_influence.csv"), row.names = FALSE)

rng <- function(x) sprintf("%.2f to %.2f (full %.2f)", min(x[-1], na.rm = TRUE),
                           max(x[-1], na.rm = TRUE), x[1])
say("\nrange across the 13 refits:")
say("  Whitsunday NTR 1987 : ", rng(inf_all$W1987))
say("  Whitsunday NTR 2004 : ", rng(inf_all$W2004))
say("  Palm NTR 1987       : ", rng(inf_all$P1987))
say("  s(maxDHW) p-value   : ", rng(inf_all$dhw_p))
say("  s(kd490) p-value    : ", rng(inf_all$kd_p))
say("  s(rugosity) p-value : ", rng(inf_all$rug_p))

# ---------------------------------------------------------------------
# 3. RESTRICTION TO THE SHELTERED STRATUM
# ---------------------------------------------------------------------
rule("3. SHELTERED SITES ONLY")
say("Adjusting for exposure assumes the functional form is right and that exposure")
say("is the only cross-sectional confounder. Restriction assumes neither. The cost")
say("is precision, and here the cost is not small.\n")
dsh <- droplevels(d[d$EXPOSURE == "Sheltered", ])
tabsh <- table(dsh$REGION[!duplicated(dsh$SITE)], dsh$NTR[!duplicated(dsh$SITE)])
say("sites in the sheltered stratum:"); cap(tabsh)
say("\nNote how thin Whitsunday's fished arm is. The restriction removes the")
say("confound but leaves few fished sites to compare against.\n")

FORM_SH <- count ~ REGION * NTR + s(YEAR, by = REGION, k = 5) +
  s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
  s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) + s(SITE, bs = "re")
mSh <- fitA(dsh, FORM_SH)
say("sheltered-only model: n = ", nrow(dsh), ", sites = ", nlevels(dsh$SITE),
    if (!is.null(mSh)) paste0(", deviance explained ",
                              sprintf("%.1f%%", summary(mSh)$dev.expl * 100)) else " (FAILED)")

# ---------------------------------------------------------------------
# 4. EACH REGION SEPARATELY
# ---------------------------------------------------------------------
rule("4. REGIONS FITTED SEPARATELY")
FORM_R <- count ~ NTR + EXPOSURE + s(YEAR, k = 5) +
  s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
  s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) + s(SITE, bs = "re")
mPalm <- fitA(droplevels(d[d$REGION == "Palm", ]), FORM_R)
mWhit <- fitA(droplevels(d[d$REGION == "Whitsunday", ]), FORM_R)
rr1 <- function(m, level) {
  if (is.null(m)) return(c(NA, NA, NA))
  b <- coef(m); V <- vcov(m); nm <- names(b); j <- match(paste0("NTR", level), nm)
  if (is.na(j)) return(c(NA, NA, NA))
  se <- sqrt(V[j, j]); unname(c(exp(b[j]), exp(b[j] - 1.96 * se), exp(b[j] + 1.96 * se)))
}
for (r in c("Palm", "Whitsunday")) {
  m <- if (r == "Palm") mPalm else mWhit
  say(r, ": n = ", sum(d$REGION == r),
      if (!is.null(m)) paste0(", deviance explained ",
                              sprintf("%.1f%%", summary(m)$dev.expl * 100)) else " (FAILED)")
}

# ---------------------------------------------------------------------
# 5. LAGGED CYCLONE, AND A DIFFERENT DISTRIBUTION
# ---------------------------------------------------------------------
rule("5. LAGGED CYCLONE AND DISTRIBUTION SWAP")
say("The Stage 1 plan asked for cyclone exposure integrated across survey intervals.")
say("That is not possible from this extract: exposure is recorded only at survey")
say("points, and the intervening years are absent. A lagged term — the previous")
say("survey's value at the same site — is what the data can support, and it is a")
say("weaker substitute. Proper integration needs external cyclone track data, which")
say("goes on the list of questions for AIMS.\n")

dlag <- droplevels(d[!is.na(d$Cyclone_lag), ])
mLag <- fitA(dlag, update(FORM, . ~ . + s(Cyclone_lag, k = 5)))
if (!is.null(mLag)) {
  say("lagged-cyclone model: n = ", nrow(dlag))
  cap(round(summary(mLag)$s.table[c("s(Cyclone)", "s(Cyclone_lag)"), , drop = FALSE], 4))
}

mTw <- tryCatch(gam(pms.leop ~ REGION * NTR + EXPOSURE + s(YEAR, by = REGION, k = 5) +
                      s(rugosity, k = 5) + s(LHC, k = 5) + s(depth, k = 5) +
                      s(kd490, k = 5) + s(maxDHW, k = 5) + s(Cyclone, k = 5) +
                      s(SITE, bs = "re"),
                    family = tw(), data = d, method = "REML"), error = function(e) NULL)
if (!is.null(mTw)) say("\nTweedie model on the density scale: deviance explained ",
                       sprintf("%.1f%%", summary(mTw)$dev.expl * 100),
                       ", p = ", sprintf("%.3f", mTw$family$getTheta(TRUE)))

# ---------------------------------------------------------------------
# 6. ALL PROTECTION ESTIMATES TOGETHER
# ---------------------------------------------------------------------
rule("6. PROTECTION UNDER EVERY SPECIFICATION")
rows <- list()
add <- function(lab, reg, lev, v) rows[[length(rows) + 1]] <<-
  data.frame(spec = lab, region = reg, protection = lev,
             ratio = v[1], lo = v[2], hi = v[3])
for (lev in c("NTR 1987", "NTR 2004")) for (reg in c("Palm", "Whitsunday")) {
  add("Full model (A)",      reg, lev, rr(mFull, reg, lev))
  add("Sheltered only",      reg, lev, rr(mSh,   reg, lev))
  add("Region fitted alone", reg, lev, rr1(if (reg == "Palm") mPalm else mWhit, lev))
  add("Lagged cyclone",      reg, lev, rr(mLag,  reg, lev))
  add("Tweedie on density",  reg, lev, rr(mTw,   reg, lev))
}
allp <- do.call(rbind, rows); rownames(allp) <- NULL
allp[, 4:6] <- round(allp[, 4:6], 2)
cap(allp[order(allp$region, allp$protection), ])
write.csv(allp, file.path(OUT, "table7_protection_all_specifications.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 7. FIGURES
# ---------------------------------------------------------------------
rule("7. FIGURES")

## Figure 10 — influence
png(file.path(OUT, "fig10_influence.png"), width = 2100, height = 1150, res = 220)
par(mfrow = c(1, 2), mar = c(3.8, 8.4, 2.6, 1.0)); base_par()
lab <- inf_all$dropped[-1]; ord <- rev(seq_along(lab))
for (cc in list(list(v = "W1987", t = "a  Whitsunday, NTR 1987"),
                list(v = "W2004", t = "b  Whitsunday, NTR 2004"))) {
  x <- inf_all[[cc$v]]; f <- x[1]; x <- x[-1]
  plot(x, ord, pch = 19, col = TEAL, cex = 1.0, yaxt = "n",
       xlim = range(c(x, f), na.rm = TRUE) * c(0.95, 1.05),
       xlab = "Density ratio vs fished", ylab = "",
       main = cc$t, font.main = 1, cex.main = 0.95, adj = 0)
  axis(2, at = ord, labels = lab, cex.axis = 0.68)
  abline(v = f, col = RED, lwd = 2)
  abline(v = 1, col = "grey70", lty = 2)
  mtext("red line = full-data estimate", side = 3, adj = 1, cex = 0.6, col = GREY)
}
dev.off(); say("wrote fig10_influence.png")

## Figure 11 — protection across all specifications
png(file.path(OUT, "fig11_protection_all_specs.png"), width = 2100, height = 1300, res = 220)
par(mfrow = c(1, 2), mar = c(3.9, 9.6, 2.6, 1.0)); base_par()
for (reg in c("Palm", "Whitsunday")) {
  s <- allp[allp$region == reg, ]
  s <- s[order(s$protection, match(s$spec, unique(allp$spec))), ]
  s$lab <- paste0(sub("NTR ", "", s$protection), " | ", s$spec)
  ok <- !is.na(s$ratio); ord <- rev(seq_len(nrow(s)))
  xr <- range(c(s$lo, s$hi, 1), na.rm = TRUE)
  plot(NA, xlim = xr, ylim = c(0.5, nrow(s) + 0.5), log = "x", yaxt = "n",
       xlab = "Density ratio vs fished", ylab = "", main = reg,
       font.main = 1, cex.main = 1.0, adj = 0)
  axis(2, at = ord, labels = s$lab, cex.axis = 0.64)
  abline(v = 1, col = "grey60", lty = 2)
  col <- ifelse(grepl("1987", s$protection), TEAL, BLUE)
  segments(s$lo[ok], ord[ok], s$hi[ok], ord[ok], col = col[ok], lwd = 2)
  points(s$ratio[ok], ord[ok], pch = 19, col = col[ok], cex = 1.0)
  if (any(!ok)) mtext("some fits did not converge", side = 1, line = 2.6,
                      adj = 1, cex = 0.6, col = RED)
}
dev.off(); say("wrote fig11_protection_all_specs.png")

writeLines(log_lines, file.path(OUT, "stage2_step3_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
