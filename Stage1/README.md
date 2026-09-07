# What explains variation in *Plectropomus leopardus* density?

Inshore reefs of the Palm and Whitsunday regions, Great Barrier Reef, 2007–2018.
An analysis of the AIMS inshore reef monitoring extract.

Ziqi (Faye) Song

---

## The result

On Whitsunday inshore reefs, sites inside no-take zones hold roughly **three times**
the density of *P. leopardus* found on fished sites. In the Palm region, 200 km away
under the same zoning framework, **no difference is detectable**.

The Whitsunday contrast survives every check applied to it. The Palm null does not
quite:

| | Full data | 13 leave-one-out refits | 5 specifications | Excludes 1? |
|---|---:|---|---|---|
| Whitsunday NTR 1987 | 3.21 | 2.88 – 3.38 | 2.58 – 3.24 | always |
| Whitsunday NTR 2004 | 2.89 | 2.59 – 3.09 | 1.87 – 2.89 | always |
| Palm NTR 1987 | 1.26 | 1.16 – 1.31 | 0.47 – 1.31 | in one of five |
| Palm NTR 2004 | 0.95 | — | 0.68 – 1.42 | never |

The five specifications are: the full model, restriction to the sheltered stratum,
each region fitted alone, a lagged cyclone term, and a Tweedie refit on the density
scale.

The one exception matters and is not buried. Restricted to sheltered sites — the
specification that handles the exposure confound by design rather than by adjustment —
Palm NTR 1987 is **0.47 (0.24–0.93)**, significantly *below* 1, where every other
specification puts it between 1.1 and 1.4. It rests on six fished and five reserve
sites. The plan's own rule for this check was to distrust the adjusted estimate when
the restriction disagrees, so the honest reading is that Palm is unresolved rather than
null. Whitsunday is unaffected: the same restriction gives 2.58 (1.56–4.25).

Nothing here is causal. Protection was not assigned at random.

## What else the analysis found

**Site identity explains more than any measured covariate.** Dropping each block from
the model in turn and recording the fall in deviance explained, the site random effect
is the largest contributor that can be separated at all. Management cannot be separated
from it — protection is fixed for a site's whole history, so removing it just hands the
work to the random effect. Measured where it is identifiable, in a model with no random
effect, management gives the largest drop of any block. Both numbers are in the report,
section 6, with bootstrap intervals. An earlier version of this analysis reported that
"management explains 48.1%"; that figure has been withdrawn as an artefact of
parameterisation — see `DECISIONS.md`.

**Habitat does not predict density.** None of rugosity (p = 0.097), live hard coral
(p = 0.375) or depth (p = 0.099) is significant, and none survives a multiplicity
correction within the habitat family. Expressed as a magnitude, rugosity briefly
appeared to exclude 1 — but splitting it into between-site and within-site components
shows neither does: 1.28 (0.94–1.75) and 1.13 (0.98–1.30). The pooled estimate looked
precise only because it averaged two weakly estimated components, and it stops
excluding 1 in any case once smoothing-parameter uncertainty is admitted.

**Something happened between 2009 and 2012 that the model cannot see.** Residuals
within a site are uncorrelated at every survey gap except one. At a three-year
separation they correlate −0.376 — a reversal, and far stronger than the −0.18 a site
random intercept induces mechanically. It is carried by the 2009-to-2012 interval in
both regions and holds under all 20 randomisations. Cyclone exposure is recorded only
at survey points, so a disturbance inside a three-year gap is invisible to the
covariate. This is that limitation made concrete rather than merely stated.

**There are more zeros than the model generates.** 52 observed against a simulated
median of 35 from the fitted model (95% interval 25–46, p = 0.004), concentrated in
Palm as the plan predicted. It is a real limitation. It does not manufacture the
protection result — excess zeros inflate apparent overdispersion, which widens
intervals rather than narrowing them, and trimming the one mostly-zero site leaves
Whitsunday at 3.27 (2.31–4.63).

**Possible residual spatial structure in Whitsunday — not established.** Moran's I on
site-level residuals is positive in Whitsunday under all 20 randomisations tested
(0.079–0.167, median 0.114) but the permutation p falls below 0.05 in only 11 of 20.
Palm is centred near zero (−0.055 to 0.026) and never significant. The residuals are
randomised quantile residuals, so anything computed from them depends on the seed; an
earlier single-seed run reported I = 0.146, p = 0.007, which overstated it. Suggestive,
not conclusive.

**Only one environmental term is left standing, and not firmly.** Intervals here use
the unconditional covariance matrix, which admits the uncertainty in the smoothing
parameters. Under it, rugosity (1.27, 0.98–1.64), turbidity (0.75, 0.55–1.01) and
cyclone exposure (0.85, 0.67–1.06) all cease to exclude 1. Only the thermal term does,
at 0.68 (0.48–0.96), and it is the only environmental term to survive a multiplicity
correction as well. Turbidity and rugosity are also unstable to which single survey is
omitted (p = 0.005–0.126 and 0.021–0.192).

**The thermal curve is only half identified.** `s(maxDHW)` is significant, but the
fitted curve rises above 4 degree heating weeks — where Whitsunday 2017 is the only
region-year with data. Refitting without that survey keeps the dip at moderate
stress and removes the rising limb, which is a spline boundary artefact. Its
worst-case concurvity against the region-specific year trends is 0.951, so it rests on
about a dozen contrasts rather than on 467 independent observations.

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
| `coral_trout_stage2_step2.R` | Habitat, environment, block contributions with bootstrap intervals. Figures 6–9. |
| `coral_trout_stage2_step3.R` | Influence analysis and sensitivity suite. Figures 10–11. |
| `coral_trout_stage2_step4.R` | Effect magnitudes, Moran's I, variable dictionary, within-between split. Figures 12–14. |
| `coral_trout_stage2_step5.R` | Concurvity, basis dimension, seed stability, depth specification, zero counts, temporal autocorrelation, multiplicity. Figures 15–16. |
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
end-to-end run from an empty `outputs/` reproduces every table, log and figure byte for
byte, verified by running it twice from empty.

Step 2 is slow — it bootstraps 200 site-level resamples at nine model fits each, which
is roughly half an hour on two cores and a few minutes on eight. Everything else takes
seconds. The resample indices are drawn up front under a fixed seed, so the answer does
not depend on how many cores do the work.

**No packages need installing.** The analysis uses base R and `mgcv`, which ships with
every R installation. That is not the same as having no versions to reconcile: mgcv's
version tracks the R version, and REML fitting and the `nb()` family have both changed
across releases, so every log records the R version, the mgcv version and the platform
that produced it. `gam()` with `s(SITE, bs = "re")` fits the random effect by REML, so
this is a mixed model in the usual sense rather than a smoothing device.

## On method

The rule followed throughout: **changes may be made in response to what is learned
about data structure or model adequacy; changes made in response to an effect
estimate are avoided.** `DECISIONS.md` logs each change against that rule, including
the cases where something was deliberately *not* changed — most notably a habitat
covariate that was not swapped after the habitat block came out null.

The project has been audited twice against its own planning documents, the second time
by an explicit hunt for claims that were stronger than the evidence. That audit
withdrew a headline number, widened every interval in the analysis, corrected two
factual errors, found a contradiction between two tables that named the same model,
and turned up three diagnostics that had never been run. All of it is in
`DECISIONS.md`, including the one audit claim that turned out to be wrong.

## Questions this analysis could not answer

- Is transect-level data with species-specific lengths available? Species-level
  biomass and size structure are the binding limitation on anything
  fisheries-relevant here.
- What survey area underlies the density values, and is the 2018 change in scaling
  convention intentional?
- Is external cyclone track data available, so exposure can be integrated across
  survey intervals rather than read at survey points? Residuals reverse sign between
  the 2009 and 2012 surveys in both regions, which says something site-differentiating
  happened in that gap that the point-sampled covariate cannot see. Track data would
  settle whether it was a storm.
- For Keppel 2021 (outside this analysis), legal density exceeds total density at 17
  of the 18 sites surveyed, by factors of 2.1 to 6.8, and implied mean weight of legal
  fish is 0.30 kg against 1.10–2.17 kg in every other region-year. Total density
  recovers cleanly on the same step, so survey area is unchanged, and the discrepancy
  is not a constant multiple. Is this known?

## Data

AIMS inshore reef monitoring, site-level extract. Not redistributed in this
repository — see `.gitignore`. Place the two source CSVs in `data/` to reproduce.
