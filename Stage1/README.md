# What explains variation in *Plectropomus leopardus* density?

Inshore reefs of the Palm and Whitsunday regions, Great Barrier Reef, 2007–2018.
An analysis of the AIMS inshore reef monitoring extract.

Ziqi (Faye) Song

---

## The result

On Whitsunday inshore reefs, sites inside no-take zones hold roughly **three times**
the density of *P. leopardus* found on fished sites. In the Palm region, 200 km away
under the same zoning framework, **no difference is detectable**.

That contrast survives every check applied to it:

| | Full data | 13 leave-one-out refits | 5 specifications | Excludes 1? |
|---|---:|---|---|---|
| Whitsunday NTR 1987 | 3.21 | 2.88 – 3.38 | 2.58 – 3.24 | always |
| Whitsunday NTR 2004 | 2.89 | 2.59 – 3.09 | 1.87 – 2.89 | always |
| Palm NTR 1987 | 1.26 | 1.16 – 1.31 | 0.47 – 1.31 | no |
| Palm NTR 2004 | 0.95 | — | 0.68 – 1.42 | no |

The five specifications are: the full model, restriction to the sheltered stratum,
each region fitted alone, a lagged cyclone term, and a Tweedie refit on the density
scale.

Nothing here is causal. Protection was not assigned at random.

## What else the analysis found

**Management explains more than site identity.** On a link-scale variance
decomposition, management takes the largest single share (0.481) against 0.089 for
unexplained site identity. Read as an upper bound — see the report, section 6.

**Habitat does not predict density.** None of rugosity (p = 0.097), live hard coral
(p = 0.375) or depth (p = 0.100) is significant. Expressed as a magnitude, rugosity
appeared to exclude 1 — but splitting it into between-site and within-site components
shows neither does: 1.28 (0.94–1.75) and 1.13 (0.98–1.30). The pooled estimate looked
precise only because it averaged two weakly estimated components.

**Possible residual spatial structure in Whitsunday — not established.** Moran's I on
site-level residuals is positive in Whitsunday under all 20 randomisations tested
(0.079–0.167, median 0.114) but the permutation p falls below 0.05 in only 11 of 20.
Palm is centred near zero (−0.055 to 0.026) and never significant. The residuals are
randomised quantile residuals, so anything computed from them depends on the seed; an
earlier single-seed run reported I = 0.146, p = 0.007, which overstated it. Suggestive,
not conclusive.

**Two environmental terms are not stable.** Turbidity ranges p = 0.005–0.126 and
rugosity p = 0.021–0.192 depending on which single survey is omitted. Neither is
reliably distinguishable from zero.

**The thermal curve is only half identified.** `s(maxDHW)` is significant, but the
fitted curve rises above 4 degree heating weeks — where Whitsunday 2017 is the only
region-year with data. Refitting without that survey keeps the dip at moderate
stress and removes the rising limb, which is a spline boundary artefact.

## Design problems the data has

**Protection is confounded with wave exposure, unevenly.** In Whitsunday, 14 of 20
fished sites are exposed while 10 of 12 long-standing reserve sites are sheltered.
In Palm the levels are comparable. The region with the large effect is the one where
the naive comparison is least trustworthy — hence the exposure adjustment and the
sheltered-only restriction.

**Thermal covariates are nearly nested in region-year.** `maxDHW` and `SSTmean` carry
99% and 98% of their variance between region-years, so they rest on about a dozen
contrasts. `Cyclone` carries 18% and is much better identified. Concurvity puts a
number on the consequence: `s(maxDHW)` has worst-case concurvity 0.951 and `s(kd490)`
0.926 against the region-specific year trends, so their coefficients are one possible
split of shared variation rather than independently identified effects.

**Depth is site-constant here**, so a penalised smooth of it was fully aliased with the
site random effect (concurvity 1.000). It is fitted as a linear term instead, which
gives identical estimates with no aliasing. Omitting it entirely raises the Whitsunday
estimate to 3.62, so it is doing real adjustment work and is kept.

**The density scaling constant is not uniform** — 0.6666 in 2007–2017, 0.6660 in 2018.
Deriving it per year recovers integers exactly. Assuming one value corrupts a year
silently.

**Species-level data is density only.** Biomass and the legal size split exist solely
pooled across all three *Plectropomus* congeners, and the mixing ratio moves between
years, so pooled biomass is not a usable proxy for one species.

## Files

| File | Contents |
|---|---|
| `coral_trout_stage1.R` | Data preparation and design audit. Figures 1–3. |
| `coral_trout_stage2_step1.R` | Negative binomial mixed model, exposure sensitivity. Figures 4–5. |
| `coral_trout_stage2_step2.R` | Habitat, environment, variance decomposition. Figures 6–9. |
| `coral_trout_stage2_step3.R` | Influence analysis and sensitivity suite. Figures 10–11. |
| `coral_trout_stage2_step4.R` | Effect magnitudes, Moran's I, variable dictionary, within-between split. Figures 12–14. |
| `coral_trout_stage2_step5.R` | Concurvity, basis-dimension check, seed stability, depth specification. Figure 15. |
| `DECISIONS.md` | Every change to the analysis plan, and what prompted it. |
| `outputs/tables.md` | All numeric tables. |
| `outputs/*_log.txt` | Full console logs — every number in the report. |
| `outputs/fig*.png` | All figures. |

## Reproducing

```bash
Rscript coral_trout_stage1.R        # must run first — writes the analysis dataset
Rscript coral_trout_stage2_step1.R
Rscript coral_trout_stage2_step2.R
Rscript coral_trout_stage2_step3.R
Rscript coral_trout_stage2_step4.R
Rscript coral_trout_stage2_step5.R
```

Run from the project directory, with the two source CSVs in `data/`. A clean
end-to-end run from an empty `outputs/` reproduces every table and log byte for byte.

**No packages are required.** The analysis uses base R and `mgcv`, both of which ship
with every R installation, so there are no versions to reconcile. `gam()` with
`s(SITE, bs = "re")` fits the random effect by REML, so this is a mixed model in the
usual sense rather than a smoothing device.

## On method

The rule followed throughout: **changes may be made in response to what is learned
about data structure or model adequacy; changes made in response to an effect
estimate are avoided.** `DECISIONS.md` logs each change against that rule, including
three cases where something was deliberately *not* changed — most notably a habitat
covariate that was not swapped after the habitat block came out null.

## Questions this analysis could not answer

- Is transect-level data with species-specific lengths available? Species-level
  biomass and size structure are the binding limitation on anything
  fisheries-relevant here.
- What survey area underlies the density values, and is the 2018 change in scaling
  convention intentional?
- Is external cyclone track data available, so exposure can be integrated across
  survey intervals rather than read at survey points?
- For Keppel 2021 (outside this analysis), legal density exceeds total density in 18
  rows and implied mean weight of legal fish is 0.30 kg against 1.11–2.17 kg
  elsewhere. Total density recovers cleanly on the same step, so survey area is
  unchanged, and the discrepancy is not a constant multiple. Is this known?

## Data

AIMS inshore reef monitoring, site-level extract. Not redistributed in this
repository — see `.gitignore`. Place the two source CSVs in `data/` to reproduce.
