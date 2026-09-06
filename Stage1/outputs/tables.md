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
observations. Cyclone exposure is not, and is much better identified.

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

## Table 4 — Protection effect across four specifications (Step 2)

Density ratio against fished sites in the same region, 95% Wald intervals.
`base` = Step 1; `+hab` adds habitat; `+env` is Model A; `sat` is Model B.

| Region | Protection | base | +hab | +env (A) | sat (B) |
|---|---|---|---|---|---|
| Palm | NTR 1987 | 1.17 | 1.29 | 1.26 | 1.27 |
| Palm | NTR 2004 | 0.96 | 1.18 | 0.95 | 0.95 |
| Whitsunday | NTR 1987 | 3.59 | 3.12 | 3.21 | 3.23 |
| Whitsunday | NTR 2004 | 3.53 | 2.81 | 2.89 | 2.90 |

Every Palm interval spans 1 in every specification. Every Whitsunday interval
excludes 1 in every specification. Adding habitat lowers the Whitsunday estimates
by 10–20%, after which they are stable.

## Table 5 — Variance decomposition on the link scale (Model A)

Unique share of linear-predictor variance.

| Block | Unique share |
|---|---:|
| Management | 0.481 |
| Space/time | 0.121 |
| Environment | 0.109 |
| Site (random effect) | 0.089 |
| Habitat | 0.047 |

Unique shares sum to 0.847; the remaining 0.153 is shared between blocks. The
largest shared component is Environment with Space/time (0.038), which is the
expected consequence of thermal covariates and year both varying at region-year
level.

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

| Region | Protection | Full (A) | Sheltered only | Region alone | Lagged cyclone | Tweedie |
|---|---|---|---|---|---|---|
| Palm | NTR 1987 | 1.26 (0.88–1.80) | **0.47 (0.24–0.93)** | 1.12 (0.75–1.67) | 1.31 (0.88–1.96) | 1.25 (0.88–1.77) |
| Palm | NTR 2004 | 0.95 (0.54–1.67) | 1.18 (0.65–2.12) | 1.42 (0.75–2.69) | 0.68 (0.37–1.25) | 1.00 (0.57–1.74) |
| Whitsunday | NTR 1987 | 3.21 (2.28–4.52) | 2.58 (1.56–4.25) | 2.73 (1.84–4.05) | 2.98 (2.10–4.25) | 3.24 (2.33–4.50) |
| Whitsunday | NTR 2004 | 2.89 (1.99–4.19) | 1.87 (1.11–3.16) | 2.49 (1.63–3.83) | 2.70 (1.84–3.96) | 2.86 (2.00–4.08) |

Every Whitsunday estimate excludes 1 in every specification. The Palm sheltered-only
estimate for NTR 1987 is the one reversal in the whole set and rests on 6 fished and
5 reserve sites; it is reported, not explained.

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

| Covariate | 10th | 90th | Ratio (95% CI) | Change |
|---|---:|---:|---|---:|
| Rugosity index | 2.44 | 4.34 | 1.27 (1.04–1.54) | +27% |
| Live hard coral (%) | 10.2 | 50.0 | 1.03 (0.81–1.32) | +3% |
| Depth (m) | 3.94 | 6.97 | 1.33 (0.95–1.88) | +33% |
| kd490 (turbidity) | 0.057 | 0.084 | 0.75 (0.57–0.97) | −25% |
| Max degree heating weeks | 0.00 | 3.38 | 0.68 (0.53–0.87) | −32% |
| Cyclone exposure index | 0.00 | 2.54 | 0.85 (0.72–0.99) | −16% |
| *Whitsunday NTR 1987 vs fished* | — | — | 3.21 (2.28–4.52) | +221% |
| *Whitsunday NTR 2004 vs fished* | — | — | 2.89 (1.99–4.19) | +189% |

The thermal contrast spans 0–3.4 DHW, the identified part of the curve, avoiding the
boundary artefact at the high end.

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

The pooled coefficient conflates two claims. Q4 asks the within-site one.

| Term | Ratio (95% CI) |
|---|---|
| Rugosity, between sites | 1.28 (0.94–1.75) |
| **Rugosity, within site** | **1.13 (0.98–1.30)** |
| Live hard coral, between sites | 1.01 (0.71–1.44) |
| Live hard coral, within site | 1.01 (0.84–1.21) |

Neither rugosity component is distinguishable from 1. The pooled rugosity estimate in
Table 8 excluded 1 only because it averaged two weakly-estimated components. **Q4's
answer is "neither"** — not structural complexity, not live coral cover.
