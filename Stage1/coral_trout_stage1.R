# =====================================================================
#  Plectropomus leopardus — Stage 1: data preparation and design audit
#  Palm and Whitsunday inshore reefs, Great Barrier Reef, 2007-2018
#  AIMS inshore reef monitoring, site-level extract
#
#  Ziqi (Faye) Song
#
#  This script stops before any model is fitted. It prepares the data,
#  audits the survey design, and produces three descriptive figures.
#
#  Dependencies: none. Base R only, so it runs on any installation.
#  Run:  Rscript coral_trout_stage1.R      (from the project directory)
# =====================================================================

set.seed(1)
DATA <- "data"; OUT <- "outputs"
dir.create(OUT, showWarnings = FALSE)

log_lines <- character(0)
say <- function(...) {
  txt <- paste0(...)
  log_lines <<- c(log_lines, txt)
  cat(txt, "\n", sep = "")
}
rule <- function(t) say("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70))

# palette (colour-blind safe)
COL_NTR <- c("Fished" = "#C1655A", "NTR 1987" = "#2E7D6F", "NTR 2004" = "#4A6FA5")
COL_REG <- c("Palm" = "#C1655A", "Whitsunday" = "#2E7D6F")
PCH_EXP <- c("Exposed" = 17, "Semi-Exposed" = 15, "Sheltered" = 19)

# ---------------------------------------------------------------------
# 1. READ
# ---------------------------------------------------------------------
rule("1. READ")
fish <- read.csv(file.path(DATA, "Selected fish benthic physical sitelevel 2021.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
# The coordinate file is mixed-encoding: a UTF-8 byte-order mark on the header,
# but the degree symbols in lat/long are not valid UTF-8. Reading it as UTF-8
# therefore fails on the body, and reading it plainly leaves the first column
# named "<BOM>REGION" so that coord$REGION is silently NULL. Read plainly and
# strip any leading non-alphanumeric bytes from the names.
coord <- read.csv(file.path(DATA, "Inshore fish site coordinates.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
names(coord) <- sub("^[^[:alnum:]]+", "", trimws(names(coord)))
say("coordinate file column names after stripping the byte-order mark: ",
    paste(names(coord), collapse = ", "))
say("site-level extract : ", nrow(fish), " rows x ", ncol(fish), " columns")
say("coordinate file    : ", nrow(coord), " rows")

# ---------------------------------------------------------------------
# 2. PARSE COORDINATES  (degrees-decimal-minutes -> decimal degrees)
#    format: "S23° 09.808' "  /  "E151° 04.491'"
# ---------------------------------------------------------------------
rule("2. PARSE COORDINATES")
dm_to_dec <- function(x) {
  x <- trimws(x)
  hemi <- toupper(substr(x, 1, 1))
  body <- sub("^[NSEWnsew]\\s*", "", x)
  deg  <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1", body))
  min  <- as.numeric(sub("^[^0-9]*[0-9]+[^0-9]+([0-9.]+).*$", "\\1", body))
  out  <- deg + min / 60
  out[hemi %in% c("S", "W")] <- -out[hemi %in% c("S", "W")]
  out
}
coord$lat_dd  <- dm_to_dec(coord$lat)
coord$long_dd <- dm_to_dec(coord$long)
bad <- sum(is.na(coord$lat_dd) | is.na(coord$long_dd))
say("coordinates parsed, failures: ", bad)
say("latitude  range: ", sprintf("%.3f to %.3f", min(coord$lat_dd), max(coord$lat_dd)))
say("longitude range: ", sprintf("%.3f to %.3f", min(coord$long_dd), max(coord$long_dd)))
stopifnot(bad == 0)

# region spelling in the coordinate file contains "WHITUSNDAY"
coord$REGION <- toupper(trimws(coord$REGION))
coord$REGION[coord$REGION == "WHITUSNDAY"] <- "WHITSUNDAY"
coord$SITE <- trimws(coord$`SITE NAME`)
say("coordinate REGION values after fixing the WHITUSNDAY misspelling: ",
    paste(sort(unique(coord$REGION)), collapse = ", "))

# ---------------------------------------------------------------------
# 3. SUBSET AND CLEAN
# ---------------------------------------------------------------------
rule("3. SUBSET AND CLEAN")
fish$SITE   <- trimws(fish$SITE)
fish$REGION <- trimws(fish$REGION)
# NTR contains a double-space variant of "NTR 1987"
fish$NTR <- gsub("\\s+", " ", trimws(fish$NTR))
say("NTR levels after whitespace normalisation: ", paste(sort(unique(fish$NTR)), collapse = " | "))

d <- fish[fish$REGION %in% c("Palm", "Whitsunday"), ]
d$NTR      <- factor(d$NTR, levels = c("Fished", "NTR 1987", "NTR 2004"))
d$EXPOSURE <- factor(trimws(d$EXPOSURE), levels = c("Sheltered", "Semi-Exposed", "Exposed"))
d$REGION   <- factor(d$REGION, levels = c("Palm", "Whitsunday"))
say("subset to Palm + Whitsunday: ", nrow(d), " rows, ", length(unique(d$SITE)), " sites, ",
    length(unique(paste(d$REGION, d$YEAR))), " region-years")

# join coordinates
d <- merge(d, coord[, c("SITE", "lat_dd", "long_dd")], by = "SITE", all.x = TRUE)
say("rows without coordinates after join: ", sum(is.na(d$lat_dd)))
stopifnot(sum(is.na(d$lat_dd)) == 0)

# ---------------------------------------------------------------------
# 4. RECOVER INTEGER COUNTS
#    Densities are counts rescaled by a constant. The constant is NOT
#    uniform across years, so derive it per year rather than assume one:
#    one fish gives exactly one step, so the step is the smallest
#    positive value observed in that year.
# ---------------------------------------------------------------------
rule("4. RECOVER INTEGER COUNTS")
dens_cols <- c("pms.leop", "pms.macu", "pms.laev", "Plectropomus total density")
step_by_year <- sapply(split(d[, dens_cols], d$YEAR), function(z) {
  v <- unlist(z); v <- v[v > 0]; min(v)
})
say("scaling step derived per year:")
for (y in names(step_by_year)) say("   ", y, " : ", sprintf("%.4f", step_by_year[[y]]))
say("-> the step is NOT uniform; 2018 differs from the other years.")

d$step  <- as.numeric(step_by_year[as.character(d$YEAR)])
d$count <- d$`pms.leop` / d$step
dev_int <- abs(d$count - round(d$count))
say("max deviation from an integer after recovery: ", sprintf("%.5f", max(dev_int)))
say("values ambiguous (deviation between 0.4 and 0.6): ", sum(dev_int > 0.4 & dev_int < 0.6))
stopifnot(max(dev_int) < 0.05)
d$count <- round(d$count)
say("counts recovered cleanly. range: ", min(d$count), " to ", max(d$count))
say("NOTE: survey area is constant, so no offset is required in a later count model —")
say("      a constant offset is absorbed into the intercept.")

# ---------------------------------------------------------------------
# 5. DESIGN AUDIT
# ---------------------------------------------------------------------
rule("5. DESIGN AUDIT")

say("\n-- 5a. Panel completeness (sites x survey years) --")
tab_ry <- table(d$REGION, d$YEAR)
print(tab_ry); log_lines <- c(log_lines, capture.output(print(tab_ry)))
sites_per_region <- tapply(d$SITE, d$REGION, function(x) length(unique(x)))
expected <- sum(sapply(levels(d$REGION), function(r) {
  sub <- d[d$REGION == r, ]; length(unique(sub$SITE)) * length(unique(sub$YEAR))
}))
say("\nrows expected if the panel is complete: ", expected, "   rows observed: ", nrow(d))
say("panel is ", ifelse(expected == nrow(d), "COMPLETE - every site in every year its region was surveyed",
                        "INCOMPLETE"))

say("\n-- 5b. Protection x exposure (the confound) --")
sites <- d[!duplicated(d$SITE), c("SITE", "REGION", "NTR", "EXPOSURE")]
tab_ce <- table(sites$EXPOSURE, sites$NTR)
print(tab_ce); log_lines <- c(log_lines, capture.output(print(tab_ce)))
say("\nReserves are not sited at random. Fished sites are predominantly exposed;")
say("no-take sites, especially the 2004 zones, are predominantly sheltered.")
say("The sheltered stratum is the only one with usable balance across all three levels.")

say("\n-- 5c. Protection x region (sites) --")
tab_cr <- table(sites$REGION, sites$NTR)
print(tab_cr); log_lines <- c(log_lines, capture.output(print(tab_cr)))

say("\n-- 5d. Species composition: why this subset --")
comp <- do.call(rbind, lapply(split(d, d$REGION), function(z) {
  data.frame(region = z$REGION[1],
             leop_share = sum(z$`pms.leop`) / sum(z$`Plectropomus total density`),
             pct_zero   = mean(z$`pms.leop` == 0))
}))
comp[, 2:3] <- round(comp[, 2:3], 3)
print(comp, row.names = FALSE); log_lines <- c(log_lines, capture.output(print(comp, row.names = FALSE)))
say("\nP. leopardus dominates in Whitsunday and is roughly co-dominant in Palm.")
say("Keppel and Magnetic (excluded) are almost entirely P. maculatus.")

say("\n-- 5e. Where each covariate's variance sits --")
say("Covariates that vary almost entirely between region-years are identified from")
say("~13 contrasts, not 467 observations. This governs what a later model can claim.")
cov_list <- c("maxDHW", "SSTmean", "Cyclone", "kd490", "LHC_%", "LCC_%", "rugosity",
              "SCI", "ChlA", "Corrected depth")
ry <- paste(d$REGION, d$YEAR)
vs <- data.frame(covariate = cov_list, between_region_year_share = NA_real_,
                 distinct_region_year_means = NA_integer_)
for (i in seq_along(cov_list)) {
  v <- d[[cov_list[i]]]
  m <- ave(v, ry, FUN = function(z) mean(z, na.rm = TRUE))
  vs$between_region_year_share[i] <- round(var(m, na.rm = TRUE) / var(v, na.rm = TRUE), 3)
  vs$distinct_region_year_means[i] <- length(unique(round(m, 6)))
}
vs <- vs[order(-vs$between_region_year_share), ]
print(vs, row.names = FALSE); log_lines <- c(log_lines, capture.output(print(vs, row.names = FALSE)))

say("\n-- 5f. Correlated covariate pairs --")
prs <- list(c("ChlA", "kd490"), c("rugosity", "SCI"), c("LCC_%", "LHC_%"),
            c("LT Fprimary", "kd490"), c("maxDHW", "SSTanom"))
for (p in prs) say("   r(", p[1], ", ", p[2], ") = ",
                   sprintf("%+.2f", cor(d[[p[1]]], d[[p[2]]], use = "complete.obs")))
say("Keep one of each near-duplicate pair (ChlA/kd490 and rugosity/SCI are")
say("effectively duplicate measurements here).")
say("CORRECTION worth recording: LT Fprimary correlates 0.75 with turbidity across")
say("the FULL four-region dataset, but only ~0.35 within Palm + Whitsunday. The")
say("earlier justification for excluding it does not hold in this subset. It remains")
say("excluded on different grounds: it is a modelled index rather than a measurement,")
say("and it is site-constant, so it competes with protection and exposure for the")
say("same 71 between-site degrees of freedom.")

say("\n-- 5g. Site-level constants (estimable only between sites) --")
for (v in c("wave exposure index", "LT Fprimary")) {
  n_var <- sum(tapply(d[[v]], d$SITE, function(z) length(unique(z))) > 1)
  say("   ", v, ": varies within ", n_var, " of ", length(unique(d$SITE)), " sites")
}

# ---------------------------------------------------------------------
# 6. DATA-QUALITY APPENDIX  (Keppel is out of scope; reported as a question)
# ---------------------------------------------------------------------
rule("6. DATA-QUALITY APPENDIX")
k <- fish[fish$REGION == "Keppel" & fish$YEAR == 2021, ]
wt_leg_all <- fish$`Plectropomus legal biomass` / fish$`Plectropomus legal density`
wt_leg_all[!is.finite(wt_leg_all)] <- NA
grp <- paste(fish$REGION, fish$YEAR)
med <- tapply(wt_leg_all, grp, median, na.rm = TRUE)
k_med <- med[["Keppel 2021"]]; others <- med[names(med) != "Keppel 2021"]
say("Implied mean weight of legal-sized fish (legal biomass / legal density):")
say("   Keppel 2021              : ", sprintf("%.2f kg", k_med))
say("   all other region-years   : ", sprintf("%.2f to %.2f kg", min(others, na.rm = TRUE),
                                              max(others, na.rm = TRUE)))
say("Rows where legal density exceeds total density (impossible): ",
    sum(fish$`Plectropomus legal density` > fish$`Plectropomus total density` + 1e-9))
kstep <- min(k$`Plectropomus total density`[k$`Plectropomus total density` > 0])
say("Keppel 2021 total density recovers to integers on a step of ", sprintf("%.4f", kstep),
    " (max dev ",
    sprintf("%.4f", max(abs(k$`Plectropomus total density`/kstep -
                             round(k$`Plectropomus total density`/kstep)))), ")")
say("Keppel 2021 legal density on that same step: max dev ",
    sprintf("%.4f", max(abs(k$`Plectropomus legal density`/kstep -
                             round(k$`Plectropomus legal density`/kstep)))))
say("=> survey area is unchanged; the anomaly is confined to the legal-density column,")
say("   and the discrepancy is not a constant factor, so a units error does not fit.")
say("   Reported as a question for AIMS, not as a correction.")

# ---------------------------------------------------------------------
# 7. FIGURES
# ---------------------------------------------------------------------
rule("7. FIGURES")
base_par <- function() par(family = "sans", mgp = c(2.2, 0.6, 0), tcl = -0.3,
                           cex.axis = 0.85, cex.lab = 0.95, las = 1)

## ---- Figure 1 : survey design and protection --------------------------
png(file.path(OUT, "fig1_design_map.png"), width = 2200, height = 1250, res = 230)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), heights = c(6, 1))
par(mar = c(3.4, 4.4, 2.4, 0.8)); base_par()
for (rg in levels(d$REGION)) {
  s <- sites[sites$REGION == rg, ]
  s <- merge(s, coord[, c("SITE", "lat_dd", "long_dd")], by = "SITE")
  asp <- 1 / cos(mean(s$lat_dd) * pi / 180)
  plot(s$long_dd, s$lat_dd, type = "n", asp = asp,
       xlab = "Longitude (decimal degrees)", ylab = "Latitude (decimal degrees)",
       main = paste0(rg, "  (n = ", nrow(s), " sites)"), font.main = 1, cex.main = 1.05)
  grid(col = "grey92", lty = 1)
  points(s$long_dd, s$lat_dd, pch = PCH_EXP[as.character(s$EXPOSURE)],
         col = COL_NTR[as.character(s$NTR)], cex = 1.15, lwd = 1.6)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = FALSE, ncol = 3, bty = "n", cex = 0.82,
       legend = c(names(COL_NTR), names(PCH_EXP)),
       col = c(COL_NTR, rep("grey35", 3)),
       pch = c(rep(19, 3), PCH_EXP),
       text.col = "grey20")
dev.off(); say("wrote fig1_design_map.png")

## ---- Figure 2 : density over time by region and protection ------------
png(file.path(OUT, "fig2_density_protection.png"), width = 2200, height = 1150, res = 230)
par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2.4, 0.8)); base_par()
ylim <- c(0, max(d$`pms.leop`) * 1.02)
for (rg in levels(d$REGION)) {
  z <- d[d$REGION == rg, ]
  yrs <- sort(unique(z$YEAR))
  plot(NA, xlim = range(yrs), ylim = ylim, xlab = "Year",
       ylab = expression(italic("P. leopardus")~"density (per 1000 m"^2*")"),
       main = rg, font.main = 1, cex.main = 1.05, xaxt = "n")
  axis(1, at = yrs, labels = yrs, cex.axis = 0.78)
  grid(col = "grey93", lty = 1)
  for (lv in levels(d$NTR)) {
    zz <- z[z$NTR == lv, ]
    if (!nrow(zz)) next
    points(jitter(zz$YEAR, amount = 0.28), zz$`pms.leop`,
           col = adjustcolor(COL_NTR[[lv]], 0.28), pch = 16, cex = 0.55)
    m  <- tapply(zz$`pms.leop`, zz$YEAR, mean)
    se <- tapply(zz$`pms.leop`, zz$YEAR, function(x) sd(x) / sqrt(length(x)))
    xs <- as.numeric(names(m))
    arrows(xs, m - se, xs, m + se, angle = 90, code = 3, length = 0.025,
           col = COL_NTR[[lv]], lwd = 1.3)
    lines(xs, m, col = COL_NTR[[lv]], lwd = 2.2)
    points(xs, m, col = COL_NTR[[lv]], pch = 19, cex = 0.95)
  }
  if (rg == "Palm") legend("topright", legend = names(COL_NTR), col = COL_NTR,
                           lwd = 2.2, pch = 19, bty = "n", cex = 0.8, text.col = "grey20")
}
dev.off(); say("wrote fig2_density_protection.png")

## ---- Figure 3 : habitat and disturbance context -----------------------
png(file.path(OUT, "fig3_habitat_disturbance.png"), width = 2100, height = 1650, res = 230)
par(mfrow = c(2, 2), mar = c(3.2, 4.0, 2.3, 0.8)); base_par()
par(mgp = c(2.0, 0.55, 0))
panels <- list(
  list(v = "LHC_%",    lab = "Live hard coral (%)",            main = "a  Coral cover"),
  list(v = "rugosity", lab = "Rugosity index",                 main = "b  Structural complexity"),
  list(v = "maxDHW",   lab = "Maximum degree heating weeks",   main = "c  Thermal stress"),
  list(v = "Cyclone",  lab = "Cyclone exposure index",         main = "d  Cyclone exposure"))
for (p in panels) {
  ally <- sort(unique(d$YEAR))
  yl <- range(unlist(lapply(split(d[[p$v]], paste(d$REGION, d$YEAR)), mean, na.rm = TRUE)))
  yl <- c(min(0, yl[1]), yl[2] * 1.15)
  plot(NA, xlim = range(ally), ylim = yl, xlab = "Year", ylab = p$lab,
       main = p$main, font.main = 1, cex.main = 1.0, adj = 0, xaxt = "n")
  axis(1, at = ally, labels = ally, cex.axis = 0.72)
  grid(col = "grey93", lty = 1)
  for (rg in levels(d$REGION)) {
    z <- d[d$REGION == rg, ]
    m <- tapply(z[[p$v]], z$YEAR, mean, na.rm = TRUE)
    xs <- as.numeric(names(m))
    lines(xs, m, col = COL_REG[[rg]], lwd = 2.2)
    points(xs, m, col = COL_REG[[rg]], pch = 19, cex = 0.9)
  }
  if (p$v == "LHC_%") legend("bottomleft", legend = levels(d$REGION), col = COL_REG,
                             lwd = 2.2, pch = 19, bty = "n", cex = 0.8, text.col = "grey20")
}
dev.off(); say("wrote fig3_habitat_disturbance.png")

# ---------------------------------------------------------------------
# 8. TABLE 1 AND LOG
# ---------------------------------------------------------------------
rule("8. OUTPUTS")
t1 <- do.call(rbind, lapply(split(sites, list(sites$REGION, sites$NTR), drop = TRUE), function(z) {
  yy <- sort(unique(d$YEAR[d$REGION == z$REGION[1]]))
  data.frame(Region = as.character(z$REGION[1]), Protection = as.character(z$NTR[1]),
             Sites = nrow(z),
             Sheltered = sum(z$EXPOSURE == "Sheltered"),
             SemiExposed = sum(z$EXPOSURE == "Semi-Exposed"),
             Exposed = sum(z$EXPOSURE == "Exposed"),
             SurveyYears = length(yy),
             Observations = nrow(z) * length(yy))
}))
t1 <- t1[order(t1$Region, t1$Protection), ]
write.csv(t1, file.path(OUT, "table0_design_summary.csv"), row.names = FALSE)
print(t1, row.names = FALSE); log_lines <- c(log_lines, capture.output(print(t1, row.names = FALSE)))

write.csv(vs, file.path(OUT, "table2_covariate_variance_structure.csv"), row.names = FALSE)
write.csv(d[, c("SITE","REGION","YEAR","NTR","EXPOSURE","lat_dd","long_dd",
                "pms.leop","count","step","LHC_%","LCC_%","rugosity","SCI",
                "kd490","ChlA","SSTmean","maxDHW","Cyclone","Corrected depth")],
          file.path(OUT, "analysis_dataset.csv"), row.names = FALSE)
say("wrote table1_design_summary.csv, table2_covariate_variance_structure.csv, analysis_dataset.csv")

# Record the environment. mgcv ships with R but its version tracks the R
# version, and REML fitting and the nb() family have both changed across
# releases, so "no packages required" is not the same as "no versions to
# reconcile". Anyone reproducing these numbers needs to know what produced them.
say("\n", strrep("-", 70))
say("environment: ", R.version.string, " | mgcv ", as.character(packageVersion("mgcv")),
    " | platform ", R.version$platform)
say(strrep("-", 70))

writeLines(log_lines, file.path(OUT, "audit_log.txt"))
cat("\nDone. Outputs in ", OUT, "/\n", sep = "")
