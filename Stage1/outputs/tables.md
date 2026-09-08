# Tables

Numeric outputs are written as CSV to `outputs/` when the scripts run, but are not
tracked in this repository. They are reproduced here so the repository is readable
without running anything.

## Table 0 — Survey design

| Region | Protection | Sites | Sheltered | Semi | Exposed | Years | Obs |
|---|---|---:|---:|---:|---:|---:|---:|
| Palm | Fished | 15 | 6 | 3 | 6 | 6 | 90 |
| Palm | NTR 1987 | 12 | 5 | 1 | 6 | 6 | 72 |
| Palm | NTR 2004 | 3 | 3 | 0 | 0 | 6 | 18 |
| Whitsunday | Fished | 20 | 3 | 3 | 14 | 7 | 140 |
| Whitsunday | NTR 1987 | 12 | 10 | 0 | 2 | 7 | 84 |
| Whitsunday | NTR 2004 | 9 | 7 | 1 | 1 | 7 | 63 |

Protection is confounded with wave exposure, and far more severely in Whitsunday
than in Palm.

## Table 2 — Where each covariate's variance sits

| Covariate | Between-region-year share | Distinct region-year means |
|---|---:|---:|
| maxDHW | 0.991 | 12 |
| SSTmean | 0.983 | 13 |
| LCC_% | 0.539 | 13 |
| LHC_% | 0.408 | 13 |
| Corrected depth | 0.298 | 2 |
| ChlA | 0.248 | 13 |
| SCI | 0.247 | 13 |
| kd490 | 0.230 | 13 |
| Cyclone | 0.182 | 7 |
| rugosity | 0.175 | 13 |

Thermal covariates are identified from roughly a dozen contrasts rather than 467
observations. Cyclone exposure is not — it is largely site-level, so better positioned than the thermal terms, which is a statement about where its variance sits rather than a guarantee it is well estimated.

## Table 3 — Protection effect, Stage 2 Step 1

Density ratio against fished sites in the same region, from a negative binomial
mixed model with a site random intercept. 95% Wald intervals.

| Region | Protection | Adjusted for exposure | Unadjusted |
|---|---|---|---|
| Palm | NTR 1987 | 1.17 (0.83–1.64) | 1.20 (0.86–1.69) |
| Palm | NTR 2004 | 0.96 (0.54–1.72) | 0.91 (0.51–1.60) |
| Whitsunday | NTR 1987 | 3.59 (2.56–5.04) | 3.29 (2.44–4.44) |
| Whitsunday | NTR 2004 | 3.53 (2.46–5.05) | 3.17 (2.29–4.41) |

Adjusting for wave exposure does not explain the regional difference. If anything
the Whitsunday estimates increase slightly once exposure is accounted for.

## Table 4 — Protection effect across five specifications (Step 2)

Density ratio against fished sites in the same region, 95% intervals on the
unconditional covariance matrix. The first rung is the Step 1 model exactly as
Step 1 fits it, so it reproduces table 3 rather than merely resembling it; depth
then enters as its own rung.

| Region | Protection | Step 1 | + depth | + hab | + env (A) | sat (B) |
|---|---|---|---|---|---|---|
| Palm | NTR 1987 | 1.17 | 1.20 | 1.29 | 1.26 | 1.27 |
| Palm | NTR 2004 | 0.96 | 0.97 | 1.18 | 0.95 | 0.95 |
| Whitsunday | NTR 1987 | 3.59 | 3.16 | 3.12 | 3.21 | 3.23 |
| Whitsunday | NTR 2004 | 3.53 | 2.95 | 2.81 | 2.89 | 2.90 |

Every Palm interval spans 1 in every specification. Every Whitsunday interval
excludes 1 in every specification.

Depth is the single covariate that moves the estimate most — Whitsunday NTR 1987
falls from 3.59 to 3.16 when it enters. An earlier version of this table folded
depth into a rung labelled "baseline (Step 1)", which put two different numbers,
3.59 and 3.16, under the same name in tables 3 and 4. Separating the rungs fixes
that and makes depth's adjustment visible.

## Table 5 — Block contributions at fixed dispersion (Step 2)

**Two earlier answers to Q2 have been withdrawn.** A variance decomposition reporting
management at 0.481, which changed to 0.391 under sum-to-zero contrasts with identical
fitted values and log-likelihood; and a deviance drop of 14.09 pp measured in a model
with no site random effect, which conflated management with everything else fixed about
a site. See `DECISIONS.md`.

Each block is dropped, the model refitted, and the fall in deviance explained recorded.
All models are fitted at the full model's theta using `negbin()` — `nb()` re-estimates
theta per model, which would put the two deviances on different scales and make the
difference uninterpretable. 95% intervals from 200 site-level cluster bootstrap
resamples. Deviance explained at fixed theta, full model: 59.5%.

| Block | Drop (pp) | 95% interval | Excludes 0 |
|---|---:|---|---|
| Site random effect (H4) | 9.76 | 1.30 – 12.43 | yes |
| Environment (H3) | 2.90 | 1.03 – 6.51 | yes |
| Regional year trends | 2.59 | 1.13 – 5.79 | yes |
| Wave exposure (confounder) | 0.23 | −1.05 – 0.94 | no |
| Habitat (H2) | −0.26 | −0.82 – 2.76 | no |
| Management (H1/Q3) | −1.69 | −3.81 – −0.44 | yes |

**Descriptive and non-additive.** The drops do not sum to anything and are not unique
variance shares.

The management row sits at or below zero because protection changes within none of the
71 sites: it and the site random intercept draw on the same between-site degrees of
freedom, so removing management hands its work to the random effect. That is not
evidence of no effect — table 7 shows the Whitsunday estimate stable across every
specification. It means deviance explained cannot separate the two, and there is no
unique, assumption-invariant share to report.

H4 expected unexplained site-level variation to be large relative to the measured
covariates. **The site random effect has the largest point estimate of any block here,
at 9.76 percentage points**, and that is the correct way to state it: its bootstrap
interval (1.30–12.43) overlaps those of environment (1.03–6.51) and the regional year
trends (1.13–5.79), so the ordering among these blocks is not established, only its
point estimate being highest. H4 is supported in direction rather than demonstrated as
a ranking.

## Table 15 — Out-of-sample prediction to unseen sites (Step 2)

**Exploratory, observational, and not causal.** What is scored here is the **marginal
predictive density of individual site-year counts at sites the model has not seen** —
marginal in the sense that the unknown site effect is integrated out rather than
estimated. It is not a variance share, it is not a test, and it carries no causal
reading.

Ten folds split by site, so a site is never in both training and test. Each model is
fitted on the training sites and scored on the held-out sites by negative binomial log
predictive density. Total log predictive density of the full model across all held-out
sites: −1214.7.

Two corrections were applied to an earlier version of this table, and both mattered:

1. **The unseen-site random effect is integrated, not set to zero.** For a site the
   model has never seen the intercept is unknown, so the predictive density is
   ∫ NB(y; exp(η + b), θ) · N(b; 0, σ²) db, not the density at b = 0 — under a log link
   the value at the median of that distribution is not its mean. The integral is taken
   by 20-node Gauss–Hermite quadrature with log-sum-exp averaging, which keeps the
   tail contributions from underflowing. σ is re-estimated inside each training fold
   from that fold's own fit, so nothing leaks from the held-out sites; across the ten
   folds it ranged 0.271 to 0.339.
2. **The standard error is clustered on sites.** The difference is summed within each
   of the 71 sites first and the spread taken across those 71 site totals. The naive
   version, which treats all 467 held-out observations as independent, is shown
   alongside because an earlier version quoted it as though it were the standard error.

| Block removed | Δ log pred. density | SE (site-clustered) | SE (naive) | Per obs | **z (clustered)** | z (naive) |
|---|---:|---:|---:|---:|---:|---:|
| Management (H1/Q3) | 48.4 | 11.7 | 10.1 | 0.1037 | **4.14** | 4.79 |
| Environment (H3) | 4.8 | 7.2 | 6.9 | 0.0103 | **0.67** | 0.7 |
| Habitat (H2) | 3 | 6.1 | 5 | 0.0065 | **0.49** | 0.6 |
| Regional year trends | 2.4 | 6.3 | 6.6 | 0.0051 | **0.38** | 0.36 |
| Wave exposure (confounder) | -7 | 3.1 | 2.1 | -0.0151 | **-2.26** | -3.33 |

The integration was much the larger of the two corrections: management's Δlpd fell from
87.3 to 48.4 and the total held-out log predictive density rose from −1232.3 to −1214.7.
Clustering the standard error contributed the smaller part, moving management's z from
4.79 to 4.14.

Knowing an unseen site's zoning improves the marginal prediction of its counts, by a
margin the other blocks do not reach in this comparison. Whether those others contribute
little or the design simply lacks the resolution to detect it is not something ten folds
over 71 sites can settle, so the right reading is that only management shows a clear
signal here — not that the rest have been shown to contribute nothing. Wave exposure
makes prediction slightly worse, which is what an adjustment term carried for confounding
control rather than for prediction can do. The two blocks that changed sign under the
correction, environment and the regional year trends, are not distinguishable from zero
either way.

**The z column is descriptive, not a formal test.** It is a ratio of a cross-validated
difference to a clustered standard error, computed on folds fixed by one seed, with the
five blocks compared and no adjustment for that. Read it as an indication of size
relative to site-to-site variability, not as a p-value in disguise.

**How not to read this.** It does not partition variance, the entries do not sum to the
model's performance, and it does not license a claim that management contributes more
variance than site identity — that would require a unique, assumption-invariant share,
which this design does not yield. It is not causal. A zoning label may predict an unseen site well because of what
protection does, or because of whatever the zoning process selected for, and this
analysis cannot tell those apart. The site random effect is deliberately absent from the
table — it is integrated out of every prediction, so its row would measure nothing and
would invite exactly the unsupported comparison above.

## Table 6 — Influence: leave one region-year out (Step 3)

Range of each estimate across 13 refits, each omitting one region-year.

| Quantity | Full data | Range across 13 refits |
|---|---:|---|
| Whitsunday NTR 1987 (ratio) | 3.21 | 2.88 – 3.38 |
| Whitsunday NTR 2004 (ratio) | 2.89 | 2.59 – 3.09 |
| Palm NTR 1987 (ratio) | 1.26 | 1.16 – 1.31 |
| s(maxDHW) p-value | <0.001 | <0.001 – <0.001 |
| s(kd490) p-value | 0.014 | 0.005 – 0.126 |
| s(rugosity) p-value | 0.097 | 0.021 – 0.192 |

The protection estimates barely move. The turbidity and rugosity terms do: both
cross conventional significance depending on which single survey is omitted, so
neither is reliably distinguishable from zero.

## Table 7 — Protection under every specification (Step 3)

Density ratio against fished sites, 95% intervals.

Four of these five were named in the project plan before anything was run. The
**lagged cyclone column is exploratory** — it was substituted after the fact for the
cyclone integration the plan asked for, which this extract cannot support because
exposure is recorded only at survey points. The plan's remaining prespecified check,
with and without Whitsunday 2017, is covered by table 6.

| Region | Protection | Full (A) | Sheltered only | Region alone | Lagged cyclone (exploratory) | Tweedie |
|---|---|---|---|---|---|---|
| Palm | NTR 1987 | 1.26 (0.87–1.82) | **0.47 (0.23–0.97)** | 1.12 (0.72–1.73) | 1.32 (0.87–2.00) | 1.25 (0.87–1.79) |
| Palm | NTR 2004 | 0.95 (0.54–1.69) | 1.18 (0.63–2.18) | 1.42 (0.72–2.80) | 0.68 (0.36–1.27) | 1.00 (0.57–1.76) |
| Whitsunday | NTR 1987 | 3.21 (2.25–4.58) | 2.58 (1.54–4.32) | 2.73 (1.83–4.08) | 2.96 (2.05–4.28) | 3.24 (2.31–4.54) |
| Whitsunday | NTR 2004 | 2.89 (1.97–4.23) | 1.87 (1.09–3.21) | 2.49 (1.61–3.86) | 2.68 (1.81–3.96) | 2.86 (1.98–4.12) |

Every Whitsunday estimate excludes 1 in every specification.

The Palm sheltered-only estimate for NTR 1987 is the one reversal in the set. It rests
on 6 fished and 5 reserve sites. Section 10 of the project plan set the rule for this
check before any of it was run: the restriction tests the same thing as the adjustment
without assuming the functional form is right, so *"if it agrees with the adjusted
estimate, say so. If it does not, distrust the adjusted one."* It does not agree.
Applying that rule, Palm is **unresolved, not null**.

Two cautions on that reversal. It is not a protection effect in the opposite direction:
eleven sites cannot distinguish a real local difference from noise, and its interval
(0.23–0.97) reaches close to 1 on the unconditional covariance matrix used throughout.
And "unresolved" is the whole claim — the analysis does not say which of the two Palm
estimates is right, only that the design cannot choose between them.

## Sheltered stratum composition

| Region | Fished | NTR 1987 | NTR 2004 |
|---|---:|---:|---:|
| Palm | 6 | 5 | 3 |
| Whitsunday | 3 | 10 | 7 |

## Table 1 — Variable dictionary (plan section 11)

Every candidate variable, its scale of variation, and why it was kept or dropped.
Selection was made in advance, not by any automatic procedure.

| Variable | Role | Varies at | Between-RY share | Status | Reason |
|---|---|---|---:|---|---|
| pms.leop / count | Response | Site-year | 0.193 | Retained | Density recovered to integer counts; area constant so no offset |
| NTR | Management | Site-constant | — | Retained | Three levels kept separate; NTR Pooled not used |
| EXPOSURE | Confounder | Site-constant | — | Retained | Confounded with protection; must accompany it |
| REGION, YEAR | Structure | — | — | Retained | Year made region-specific after Step 1 residual diagnostics |
| SITE | Structure | — | — | Retained | Random intercept, 71 levels, REML |
| rugosity | Habitat | Site-year | 0.175 | Retained | Structural complexity; chosen over SCI (r = 0.93) |
| LHC_% | Habitat | Site-year | 0.408 | Retained | Live hard coral; chosen over LCC_% (r = 0.80) in advance |
| Corrected depth | Habitat | Site-constant here | 0.298 | Retained | Varies at 0 of 71 sites in this subset |
| kd490 | Environment | Site-year | 0.230 | Retained | Water clarity; chosen over ChlA (r = 1.00 in this subset) |
| maxDHW | Disturbance | Region-year | 0.991 | Retained (Model A only) | Nearly nested in region-year; not identifiable in Model B |
| Cyclone | Disturbance | Mostly site-level | 0.182 | Retained | Better identified than maxDHW; interval integration impossible here |
| SCI | Habitat | Site-year | 0.247 | Dropped | Duplicate of rugosity (r = 0.93) |
| ChlA | Environment | Site-year | 0.248 | Dropped | Duplicate of kd490 (r = 1.00) |
| LCC_% | Habitat | Site-year | 0.539 | Dropped | Correlates 0.80 with LHC_%; **not** reinstated after the habitat null |
| SSTmean | Environment | Region-year | 0.983 | Dropped | 98% between-region-year; would compete with maxDHW |
| LT Fprimary | Fishing proxy | Site-constant | — | Dropped | Modelled index, and site-constant, so competes with protection for the same 71 df |
| wave exposure index | Confounder | Site-constant | — | Dropped | Redundant with categorical EXPOSURE |

## Table 8 — Effect magnitude across the observed range (Step 4)

Each covariate moved from its 10th to its 90th percentile, everything else held fixed.
Intervals use the **unconditional** covariance matrix, which admits the uncertainty in
the smoothing parameters; the conditional default treats them as known and is too narrow
for anything read off a smooth.

| Covariate | 10th | 90th | Ratio | 95% CI (unconditional) | Change |
|---|---:|---:|---:|---|---:|
| Rugosity index | 2.44 | 4.336 | 1.27 | 0.98–1.64 | +27% |
| Live hard coral (%) | 10.24 | 50 | 1.03 | 0.79–1.34 | +3% |
| Depth (m) | 3.94 | 6.97 | 1.33 | 0.94–1.89 | +33% |
| kd490 (turbidity) | 0.057 | 0.084 | 0.74 | 0.55–1.01 | -25% |
| Max degree heating weeks | 0 | 3.38 | 0.68 | 0.48–0.96 | -32% |
| Cyclone exposure index | 0 | 2.535 | 0.84 | 0.67–1.06 | -16% |
| *Whitsunday NTR 1987 vs fished* | — | — | 3.21 | 2.25–4.58 | +221% |
| *Whitsunday NTR 2004 vs fished* | — | — | 2.89 | 1.97–4.23 | +189% |

Rugosity, turbidity and cyclone exposure do not exclude 1. Only the thermal term does,
and it is the only environmental term to survive the multiplicity correction in table 14
as well — which makes it the most suggestive of them rather than an established effect.
It rests on roughly a dozen region-year contrasts, its worst-case concurvity against the
region-specific year trends is 0.951, and those year smooths carry a basis warning that
does not clear at k = 6 (table 18). Switching to the unconditional matrix inflates the standard error by 33–42% for
terms read off a smooth and by about 1% for the parametric protection contrasts, which is
why the headline result is untouched by the change.

The thermal contrast spans 0–3.4 DHW, the identified part of the curve, avoiding the
boundary artefact at the high end. Its worst-case concurvity against the region-specific
year trends is 0.951, so it is one possible split of shared variation rather than an
independently identified effect.

## Table 9 — Spatial autocorrelation in site-level residuals (Step 4)

Inverse-distance weights, 4999 permutations, tested within region.

Single-seed values (seed 1), **superseded by Table 11** — these depend on the
randomisation in the quantile residuals.

| Region | Sites | Moran's I | Expected | p | Median separation |
|---|---:|---:|---:|---:|---:|
| Palm | 30 | 0.011 | −0.034 | 0.878 | 5.4 km |
| Whitsunday | 41 | 0.146 | −0.025 | 0.007 | 8.9 km |

## Table 11 — Seed stability of the residual diagnostics (Step 5)

Randomised quantile residuals draw from a uniform, so every diagnostic built on them
depends on the seed. Twenty seeds:

| | Range | Median | p < 0.05 |
|---|---|---:|---|
| Whitsunday Moran's I | 0.079 – 0.167 | 0.114 | 11 of 20 |
| Palm Moran's I | −0.055 – 0.026 | −0.017 | 0 of 20 |
| Shapiro–Wilk on residuals | p 0.298 – 0.995 | — | 0 of 20 rejected |

**Correction.** Table 9's single-seed result overstated the Whitsunday finding. What is
defensible: Moran's I is positive in Whitsunday under every seed and near zero in Palm
under every seed, but the Whitsunday value is not significant at conventional thresholds
under half of them. Residual normality, by contrast, is robust.

## Table 12 — Depth specification (Step 5)

| Specification | AIC | Deviance explained | Whitsunday NTR 1987 |
|---|---:|---:|---|
| depth as a smooth | 2353.1 | 60.2% | 3.21 (2.28–4.52) |
| **depth linear (used)** | 2353.0 | 60.2% | 3.21 (2.28–4.52) |
| depth omitted | 2354.1 | 59.3% | 3.62 (2.62–4.99) |

Depth is constant within every site here, so a penalised smooth of it was fully aliased
with the site random effect (concurvity 1.000). The linear form gives identical estimates
with no aliasing. Omitting it raises the Whitsunday estimate, so it is kept.

## Concurvity (Step 5)

Worst-case concurvity against the rest of the model:

| Term | Worst |
|---|---:|
| s(SITE) | 1.000 |
| s(maxDHW) | 0.951 |
| s(kd490) | 0.926 |
| s(YEAR):Whitsunday | 0.914 |
| s(YEAR):Palm | 0.879 |
| s(LHC) | 0.750 |
| s(Cyclone) | 0.655 |
| s(rugosity) | 0.644 |

The thermal and turbidity smooths are substantially reproducible from the year trends.
Their coefficients are one possible split of shared variation, not independently
identified effects.

## Table 10 — Habitat effects split within and between sites (Step 4)

A pooled habitat coefficient conflates two claims: that sites with more structure hold
more fish, and that a site which gains structure gains fish. Splitting each covariate
into a site mean and a within-site deviation separates them. Intervals are unconditional.

| Term | 10th | 90th | Ratio (95% CI) |
|---|---:|---:|---|
| Rugosity, between sites | 2.827 | 3.873 | 1.28 (0.88–1.85) |
| Rugosity, within site | -0.797 | 0.674 | 1.13 (0.98–1.31) |
| Live hard coral, between sites | 15.733 | 44.171 | 1.01 (0.69–1.47) |
| Live hard coral, within site | -15.537 | 14.046 | 1.01 (0.81–1.24) |

Neither component of either habitat variable excludes 1. The pooled rugosity estimate
looked precise only because it averaged two weakly estimated components — and on the
unconditional covariance matrix the pooled estimate does not exclude 1 either. Q4 asks
whether the habitat association runs with structural complexity, with live coral cover,
or with neither. On these estimates, no clear association is found for either variable —
which answers Q4 for the two candidates the plan named, and is not a general finding
about habitat.

## Table 13 — Goodness-of-fit bootstrap (Step 5)

Two questions from one parametric bootstrap, 300 replicates. Each replicate draws a
fresh set of site effects from N(0, 0.3065²), simulates a dataset, **refits the same
model to it**, and recomputes the statistic from the refit exactly as it is computed
from the real fit.

Three things in that sentence are corrections to an earlier version. It refits, rather
than simulating at fixed means, so parameter estimation is inside the reference
distribution. It redraws the site effects rather than reusing the fitted ones, which had
understated between-site variation — the fitted site effects have SD 0.219 against an
estimated component of 0.307 — and made the null too narrow. And the residual-correlation
statistic is the most negative correlation across **all eleven** survey gaps, so the
search is part of the statistic rather than a choice made after seeing the answers.

| Statistic | Observed | Simulated median | 95% interval | p | B |
|---|---:|---:|---|---:|---:|
| Zero discrepancy (observed − expected zeros) | 16.45 | 1.22 | −6.58 – 9.24 | **0.0033** | 300 |
| Most negative gap correlation | −0.376 | −0.259 | −0.452 – −0.176 | 0.0797 | 300 |

p is computed as (extreme + 1) / (B + 1) and is never zero; the smallest attainable
value with 300 replicates is 0.0033.

**Zeros: a real excess.** About 16 more zeros than the model expects. The negative
binomial is not a complete description of the zero process.

On consequences, the narrow claim only: the two sensitivities in table 16 did not
materially change the estimate, and there is a general argument that a zero excess
inflates apparent overdispersion and so widens intervals rather than narrowing them.
**No hurdle or zero-inflated negative binomial was fitted**, and that refit is what
would settle it. What is supported is that the checks that were run did not move the
estimate — not that the misfit is harmless.

**Residual temporal correlation: unresolved, and exploratory.** The observed value sits
toward the extreme of the simulated null without being clearly outside it, so the check
neither confirms nor rules out residual temporal structure. It is left unresolved rather
than pushed to one side of a threshold, and the exact p is recorded in the table above
rather than carried into the narrative, because reading it as a decision would be to
treat an exploratory diagnostic as a formal test.

An earlier version reported this as a clear finding and offered an explanation for it.
That came from testing against zero rather than against the correlation a site random
intercept induces, from selecting the strongest of eleven gaps after seeing them, and
from simulating without redrawing the site effects. All three made the test more
permissive. No mechanism is proposed for what remains: this analysis has no covariate
that would distinguish one, and the survey design cannot separate candidates. More
replicates would sharpen the estimate but change no conclusion here, so none were added.

Observed correlations by gap, for the record:

| Gap (years) | 1 | 2 | **3** | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| r | 0.046 | −0.028 | **−0.376** | −0.090 | −0.085 | −0.194 | −0.128 | 0.069 | −0.028 | 0.126 | −0.087 |

## Table 16 — Does either diagnostic move the estimate? (Step 5)

| Model | Whitsunday NTR 1987 | Whitsunday NTR 2004 |
|---|---|---|
| Full model (A) | 3.21 (2.25–4.58) | 2.89 (1.97–4.23) |
| Zero-trimmed | 3.27 (2.31–4.63) | 2.89 (1.98–4.22) |
| + site-year random slope | 3.27 (2.10–5.08) | 2.90 (1.79–4.71) |

The random-slope row is the sensitivity the plan called for: sites are allowed their own
trajectory rather than only their own level. The wider interval is the cost of the extra
flexibility, not a change in the estimate.

## Table 17 — Observed against fitted, final model (Step 5)

Step 1 checked this for the baseline; the final model had never been checked the same
way. Deciles of fitted value, with the standard error of the observed mean in each bin.

The **z column is descriptive, not a test.** The bins are defined by the fitted values
themselves, the ten are not independent, and no multiplicity or selection adjustment is
applied — so |z| > 2 in one bin is not a rejection. Read it as a scale for how far each
bin sits from its prediction relative to the noise in that bin, and read the pattern
across bins rather than any single value.

| Bin (fitted) | n | Fitted mean | Observed mean | SE | z |
|---|---:|---:|---:|---:|---:|
| 0.57 – 1.80 | 47 | 1.34 | 0.91 | 0.21 | −2.05 |
| 1.80 – 2.38 | 47 | 2.11 | 1.70 | 0.24 | −1.71 |
| 2.38 – 3.06 | 46 | 2.72 | 2.72 | 0.27 | 0.00 |
| 3.06 – 3.54 | 47 | 3.27 | 3.30 | 0.40 | 0.07 |
| 3.54 – 4.32 | 47 | 3.92 | 3.91 | 0.40 | −0.02 |
| 4.32 – 5.57 | 46 | 5.00 | 5.43 | 0.39 | 1.10 |
| 5.57 – 7.04 | 47 | 6.39 | 6.87 | 0.45 | 1.07 |
| 7.04 – 9.06 | 46 | 8.02 | 7.76 | 0.63 | −0.41 |
| 9.06 – 11.8 | 47 | 10.21 | 10.64 | 0.64 | 0.67 |
| 11.8 – 24.8 | 47 | 15.64 | 15.98 | 0.86 | 0.40 |

The lowest bin is the one to watch: the model over-predicts there (z = −2.05), which is
the zero excess of table 13 seen from another angle rather than a separate problem.
Elsewhere the mean structure tracks the data closely, with no systematic drift across
the range.

## Table 18 — Basis dimension, tested rather than asserted (Step 5)

An earlier version stated that the year-smooth flag could not be fixed by raising k and
that the kd490 flag was a real limitation. Neither had been checked. Both are refitted
here at the largest basis the data allow.

| Term | Specification | k-index | p | edf | AIC | Whitsunday NTR 1987 |
|---|---|---:|---:|---:|---:|---|
| s(kd490) | k = 5 (used) | 0.847 | 0.0025 | 2.68 | 2353.0 | 3.21 (2.25–4.58) |
| s(kd490) | k = 10 | 0.847 | 0.0025 | 2.79 | 2353.2 | 3.21 (2.25–4.58) |
| s(kd490) | k = 20 | 0.847 | 0.0025 | 2.81 | 2353.2 | 3.21 (2.25–4.58) |
| s(YEAR):Palm | k = 5 (used) | 0.862 | 0.0075 | 3.46 | 2353.0 | 3.21 (2.25–4.58) |
| s(YEAR):Palm | k = 6 (maximum) | 0.863 | 0.0100 | 4.02 | 2344.8 | 3.18 (2.24–4.52) |
| s(YEAR):Whitsunday | k = 5 (used) | 0.862 | 0.0050 | 1.34 | 2353.0 | 3.21 (2.25–4.58) |
| s(YEAR):Whitsunday | k = 6 (maximum) | 0.863 | 0.0075 | 1.01 | 2344.8 | 3.18 (2.24–4.52) |

Neither flag moves with k. Raising the kd490 basis from 5 to 20 leaves the k-index at
0.847 to three decimals, and the year smooths cannot exceed k = 6 because Palm has only
six distinct survey years. So in both cases the residual pattern is not something a
larger basis can represent, and in both cases the protection estimate is unmoved. Model
B, which saturates time with region-year as a factor, is the specification that answers
the year flag directly rather than by enlarging a basis, and it agrees with Model A on
protection to within 0.02. The kd490 flag remains a real limitation on what can be
claimed about turbidity — one of several reasons turbidity is reported as unresolved.

## Table 14 — Multiplicity within hypothesis families (Step 5)

Model A reports 17 tests. One correction across all 17 would be the wrong instrument,
since they answer four different questions. Holm within each of the four hypothesis
families the plan set out in advance. Holm controls the family-wise error rate, needs no
independence assumption — which matters here, since concurvity between these terms is
high — and is uniformly at least as powerful as Bonferroni.

| Family | Term | p | p (Holm) | Survives |
|---|---|---:|---:|---|
| H1/Q3 management | NTR 1987 | 0.2099 | 0.4198 | no |
| H1/Q3 management | NTR 2004 | 0.8629 | 0.8629 | no |
| H1/Q3 management | REGION x NTR 1987 | 0.0004 | 0.0014 | **yes** |
| H1/Q3 management | REGION x NTR 2004 | 0.0007 | 0.0020 | **yes** |
| H2 habitat | s(rugosity) | 0.0971 | 0.2912 | no |
| H2 habitat | s(LHC) | 0.3748 | 0.3748 | no |
| H2 habitat | depth | 0.0991 | 0.2912 | no |
| H3 environment | s(maxDHW) | <0.0001 | <0.0001 | **yes** |
| H3 environment | s(kd490) | 0.0136 | 0.0272 | **yes** |
| H3 environment | s(Cyclone) | 0.0776 | 0.0776 | no |
| H4 site variation | s(SITE) | <0.0001 | <0.0001 | **yes** |

Q3, the region-by-protection interaction, clears the correction by an order of
magnitude. The NTR main effects do not, but those are the Palm estimates, which were
already unresolved; failing a multiplicity correction does not make an unresolved result
resolved in either direction.

Nothing in the habitat family survives, which agrees with the within-between split in
table 10.

The environment family is where two defensible corrections part company. The thermal
term survives Holm *and* the unconditional interval in table 8. Cyclone exposure
survives neither. Turbidity survives Holm but not the interval, and its basis-dimension
flag in table 18 does not clear at any k, so it is reported as unresolved — when two
corrections disagree, the weaker reading is the one quoted.

The year-trend smooths are in no family: they are adjustment terms rather than
hypotheses, and correcting a hypothesis test for the significance of a nuisance term
would be a category error.
