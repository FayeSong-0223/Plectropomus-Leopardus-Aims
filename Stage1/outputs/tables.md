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

## Table 6 — Influence: leave one region-year out (Step 3)

Range of each estimate across 13 refits, each omitting one region-year.

| Quantity | Full data | Range across 13 refits |
|---|---:|---|
| Whitsunday NTR 1987 (ratio) | 3.21 | 2.91 – 3.40 |
| Whitsunday NTR 2004 (ratio) | 2.89 | 2.63 – 3.10 |
| Palm NTR 1987 (ratio) | 1.26 | 1.16 – 1.30 |
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
