# What explains variation in *Plectropomus leopardus* density?

Inshore reefs of the Palm and Whitsunday regions, Great Barrier Reef, 2007–2018.
An analysis of the AIMS inshore reef monitoring extract.

Ziqi (Faye) Song

---

## The result

On Whitsunday inshore reefs, sites inside no-take zones hold roughly **three times**
the density of *P. leopardus* found on fished sites. In the Palm region, 200 km away
under the same zoning framework, **no difference is detectable and the specifications
disagree**.

This is a robust association, not a measured effect of protection. Every specification
tried leaves the Whitsunday estimate where it is, and adjusting for wave exposure or
restricting to sheltered sites does not weaken it — so wave exposure is not the
explanation. Reserve placement was not random, though, and no observational check here
can establish that whatever else distinguishes the zoned sites is innocuous.

Every check applied leaves the Whitsunday estimate where it is. Palm does not hold
together as a null:

| | Full data | 13 leave-one-out refits | 5 specifications | Excludes 1? |
|---|---:|---|---|---|
| Whitsunday NTR 1987 | 3.21 | 2.88 – 3.38 | 2.58 – 3.24 | always |
| Whitsunday NTR 2004 | 2.89 | 2.59 – 3.09 | 1.87 – 2.89 | always |
| Palm NTR 1987 | 1.26 | 1.16 – 1.31 | 0.47 – 1.31 | in one of five |
| Palm NTR 2004 | 0.95 | — | 0.68 – 1.42 | never |

The five specifications are: the full model, restriction to the sheltered stratum,
each region fitted alone, a lagged cyclone term, and a Tweedie refit on the density
scale. Four of those were named in the project plan before anything was run. The
**lagged cyclone term was not** — it is exploratory, substituted after the fact for the
cyclone integration the plan asked for, which this extract cannot support because
exposure is recorded only at survey points.

The one exception matters and is not buried. Restricted to sheltered sites — the
specification that handles the exposure confound by design rather than by adjustment —
Palm NTR 1987 is **0.47 (0.23–0.97)**, below 1 though its interval reaches close to
it, where every other specification puts it between 1.1 and 1.4. It rests on six fished and five reserve
sites. The plan's own rule for this check was to distrust the adjusted estimate when
the restriction disagrees, so the honest reading is that Palm is unresolved rather than
null. Whitsunday is unaffected: the same restriction gives 2.58 (1.54–4.32).

Nothing here is causal. Protection was not assigned at random.

## What else the analysis found

**How much does management explain? Not as a unique variance share.** Protection
changes within 0 of 71 sites, so its entire contribution is between-site — and a site
random intercept is also entirely between-site. The two draw on the same degrees of
freedom, and the consequence is specific: **there is no unique, assumption-invariant
variance share to report.** Any number that appears to separate them is doing so on the
strength of an assumption the data do not test — a choice of contrast coding, or the
removal of the random effect — and a different, equally defensible assumption gives a
different number.

That is a statement about what this design identifies, not a claim that the question is
unanswerable in principle: a design with protection changing within sites, or with
external information about how sites were selected for zoning, could address it. Two
earlier attempts here have been withdrawn — a variance share of 48.1% that depended on
the coding of the factor contrasts, and a deviance drop measured in a model with no
random effect, which conflated management with everything else fixed about a site. Both
are documented in `DECISIONS.md`.

What is reported instead is a drop-one-block comparison at fixed dispersion —
descriptive, non-additive, not variance shares — and an exploratory out-of-sample
comparison. On the first, the **site random effect has the largest point estimate** of
any block, at 9.76 percentage points of deviance explained; its bootstrap interval
(1.30–12.43) overlaps those of environment (1.03–6.51) and the regional year trends
(1.13–5.79), so the ordering among blocks is not established and only the point estimate
is highest. H4 is supported in direction rather than demonstrated as a ranking. That second one predicts **individual site-year counts at sites the model
has never seen, marginally over the unknown site effect**. It is exploratory,
observational and not causal, and stays that way whatever its z.

Both defects an earlier version of that comparison carried have now been corrected. The
unseen-site random effect is integrated over rather than set to zero — under a log link
the prediction at the median of the random-effect distribution is not its mean, so the
integral is the right quantity, taken by 20-node Gauss–Hermite quadrature with
log-sum-exp averaging, with σ re-estimated inside each training fold (range 0.271–0.339
across the ten). And the standard error is now clustered on the 71 sites rather than
treating all 467 held-out observations as independent.

| Block removed | Δ log pred. density | SE (site) | **z** |
|---|---:|---:|---:|
| Management | **48.4** | 11.7 | **4.14** |
| Environment | 4.8 | 7.2 | 0.67 |
| Habitat | 3.0 | 6.1 | 0.49 |
| Regional year trends | 2.4 | 6.3 | 0.38 |
| Wave exposure | −7.0 | 3.1 | −2.26 |

The integration was much the larger correction — management's Δlpd fell from 87.3 to
48.4; clustering moved its z from 4.79 to 4.14. Knowing an unseen site's zoning improves
the marginal prediction of its counts, by a margin the other blocks do not reach here;
whether those others contribute little or the design simply lacks the resolution to
detect it is not something ten folds over 71 sites can settle.

**The z is descriptive, not a formal test.** It is a ratio of a cross-validated
difference to a bootstrap-free clustered standard error, computed on folds chosen by
one seed, with the blocks selected in advance but no adjustment for comparing five of
them. Read it as an indication of how large the difference is relative to its
site-to-site variability, not as a p-value in disguise.

The result is a predictive statement and nothing more. It does not partition variance,
it does not license a claim that management contributes more variance than site
identity, and it is not causal: a zoning label may predict an unseen site well because
of what protection does, or because of whatever the zoning process selected for, and
this analysis cannot tell those apart.

**Habitat does not predict density.** None of rugosity (p = 0.097), live hard coral
(p = 0.375) or depth (p = 0.099) is significant, and none survives a multiplicity
correction within the habitat family. Expressed as a magnitude, rugosity briefly
appeared to exclude 1 — but splitting it into between-site and within-site components
shows neither does: 1.28 (0.88–1.85) and 1.13 (0.98–1.31). The pooled estimate looked
precise only because it averaged two weakly estimated components, and it stops
excluding 1 in any case once smoothing-parameter uncertainty is admitted.

**The negative binomial does not fully account for the zeros.** Against a parametric
bootstrap that refits every replicate and redraws the site effects from their estimated
distribution, the data hold about 16 more zeros than the model expects, against a
simulated distribution centred near 1 (p = 0.003, 300 replicates). It is a limitation on
the distributional form.

What can be said about its consequences is narrower than it might appear. The two
sensitivities available here — trimming the mostly-zero site, and giving sites their own
time trajectory — did not materially change the estimate, and there is a general argument
that a zero excess inflates apparent overdispersion and so widens intervals rather than
narrowing them. But **no hurdle or zero-inflated negative binomial was fitted**, and that
is the refit that would actually settle it. The honest statement is that the checks that
were run did not move the estimate, not that the misfit is harmless.

**Residual temporal correlation: unresolved, and treated as exploratory.** The most
negative residual correlation across the eleven survey gaps is −0.376, which sits inside
the range a correctly specified model produces by itself once the site random effect is
redrawn — a simulated null centred on −0.259 with a 95% interval of −0.452 to −0.176.
The observed value is toward the extreme of that null without being clearly outside it,
so the check neither confirms nor rules out residual temporal structure and is left
unresolved rather than pushed to one side.

An earlier version reported this as a clear finding and offered an explanation for it.
That came from testing the raw correlation against zero rather than against the
correlation a site random intercept induces by construction, from selecting the
strongest of eleven gaps after seeing them, and from simulating without redrawing the
site effects. Each made the test more permissive. No mechanism is proposed for what
remains: there is no covariate here that would distinguish one.

**Neither diagnostic moves the protection estimate.** Trimming the mostly-zero site gives
3.27 (2.31–4.63); adding a site-level random slope on year gives 3.27 (2.10–5.08), the
wider interval being the cost of the extra flexibility rather than a change in the
estimate.

**Possible residual spatial structure in Whitsunday — not established.** Moran's I on
site-level residuals is positive in Whitsunday under all 20 randomisations tested
(0.079–0.167, median 0.114) but the permutation p falls below 0.05 in only 11 of 20.
Palm is centred near zero (−0.055 to 0.026) and never significant. The residuals are
randomised quantile residuals, so anything computed from them depends on the seed; an
earlier single-seed run reported I = 0.146, p = 0.007, which overstated it. Suggestive,
not conclusive.

**One environmental term is suggestive; none is established.** Intervals here use the
unconditional covariance matrix, which admits the uncertainty in the smoothing
parameters. Under it, rugosity (1.27, 0.98–1.64), turbidity (0.75, 0.55–1.01) and
cyclone exposure (0.85, 0.67–1.06) all cease to exclude 1. Only the thermal term does,
at 0.68 (0.48–0.96), and it is the only environmental term to survive the multiplicity
correction as well.

**That makes it suggestive, not established**, and three things stand between it and a
firmer reading. `maxDHW` carries 99% of its variance between region-years, so the term
rests on roughly a dozen region-year contrasts rather than on 467 observations. Its
worst-case concurvity against the region-specific year trends is 0.951, so its
coefficient is one possible split of variation those terms share rather than an
independently identified effect. And those year smooths carry a basis-dimension warning
that does not clear at k = 6, the largest the survey years allow — so one of the terms
competing for that shared variation is itself unresolved. A negative association at
moderate thermal stress is consistent with these data and worth pursuing with a design
that can identify it; it is not demonstrated here.

Turbidity is unresolved on three separate grounds: its interval includes 1, its p-value
ranges 0.005–0.126 depending on which single survey is omitted, and its own basis flag
does not clear at k = 20. Rugosity is likewise unstable (p = 0.021–0.192).

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
| `coral_trout_stage2_step5.R` | Concurvity, basis dimension, seed stability, depth specification, goodness-of-fit bootstrap, multiplicity, observed vs fitted. Figures 15–16. |
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

Two steps are slow. Step 2 bootstraps 200 site-level resamples at seven model fits each
and then cross-validates over ten site-blocked folds; Step 5 simulates and refits 300
times. Together they take roughly 50 minutes on two cores and a good deal less on eight.
Everything else takes seconds. All resample indices and simulation seeds are drawn up
front under fixed seeds, so the answer does not depend on how many cores do the work.

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
  is not a constant multiple. Is this known? (The whole extract holds 18 such rows: the
  eighteenth is Whitsunday 2018 site HY3, excess 0.0018 — a rounding artefact, not the
  same thing.)

## Data

AIMS inshore reef monitoring, site-level extract. Not redistributed in this
repository — see `.gitignore`. Place the two source CSVs in `data/` to reproduce.
