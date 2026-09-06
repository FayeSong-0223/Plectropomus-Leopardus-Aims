# Tables

Numeric outputs are written as CSV to `outputs/` when the scripts run, but are not
tracked in this repository. They are reproduced here so the repository is readable
without running anything.

## Table 1 — Survey design

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
