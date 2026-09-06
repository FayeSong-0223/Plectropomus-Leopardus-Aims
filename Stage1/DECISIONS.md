# Decisions log

Every change to the analysis plan, with what prompted it. Kept so that "the final
model differs from the plan" can be answered with a reason rather than a shrug.

The rule being followed: changes may be made in response to what is learned about
**data structure or model adequacy**. Changes made in response to an **effect
estimate** are avoided, because that is how a model gets chosen for the answer it
gives.

---

**2026-09-05 — Species and region restricted before any analysis.**
*P. leopardus* only; Palm and Whitsunday only. Keppel and Magnetic are almost
entirely *P. maculatus*, so including them would confound species with region.
Prompted by data structure, decided before any model was fitted.

**2026-09-05 — Scaling constant derived per year rather than assumed.**
The constant is 0.6666 in 2007–2017 and 0.6660 in 2018. Assuming a single value
corrupts one year's counts silently. Prompted by a consistency check.

**2026-09-05 — No offset in the count model.**
Survey area is constant across all observations, so a constant offset is absorbed
into the intercept and changes no covariate coefficient.

**2026-09-05 — `LT Fprimary` excluded, but the stated reason was corrected.**
Originally excluded because it correlates 0.75 with turbidity. That figure is from
the full four-region dataset; within Palm and Whitsunday it is only about 0.35, so
the collinearity argument does not hold here. It remains excluded on different
grounds: it is a modelled index rather than a measurement, and it is constant
within every site, so it competes with protection and exposure for the same 71
between-site degrees of freedom.

**2026-09-06 — `mgcv` chosen over `glmmTMB` for the mixed model.**
`gam()` with `s(SITE, bs = "re")` and `family = nb()` fits the random effect by
REML, so it is a mixed model in the usual sense. `mgcv` ships with every R
installation, and Step 2 needs it for smooth terms anyway, so the whole analysis
stays in one model class. Prompted by tooling, not by results.

**2026-09-06 — Habitat covariates will be adjusted for in Step 2, decided in advance.**
If reserves were *placed* in different habitat, adjusting is correct. If protection
*changed* the habitat, habitat is a mediator and adjusting would remove part of the
effect being measured. For coral trout over this period the first is far more
plausible. Recorded before fitting so the choice is not revisited after seeing which
way the coefficient moves.

**2026-09-06 — Year will become region-specific in Step 2.**
Step 1 used a single `factor(YEAR)` shared across regions. Residuals by region-year
(Figure 4c) show clear structure: Palm 2014–2018 sits consistently below zero and
Palm 2007–2012 above it, while Whitsunday does not follow the same pattern. The two
regions have different temporal trajectories and a shared year term cannot represent
them. Prompted by a diagnostic, not by an effect estimate.

**2026-09-06 — Two models fitted rather than one, to resolve a genuine conflict.**
Step 1's diagnostics required region-specific year effects. But region-year fixed
effects absorb 99% of `maxDHW`'s variance, making a thermal effect inestimable.
Rather than choose silently: Model A uses region-specific *smooth* year trends,
leaving year-to-year deviation available to the environmental covariates; Model B
saturates time with region-year as a factor and drops the thermal terms as
unidentifiable. Environmental effects are read from A only, protection from both.
The two agree on protection to two decimal places.

**2026-09-06 — Thermal sensitivity check added after seeing the fitted curve.**
`s(maxDHW)` was the strongest environmental term, but the fitted curve *rose* above
about 4 degree heating weeks, which would mean heat stress increases coral trout
density. The rug showed a gap between 4 and 5.8 DHW, and Whitsunday 2017 is the only
region-year above 4 (41 of 467 observations). Model A was refitted without that
survey. This is a change prompted by an implausible fitted shape and a data-coverage
gap — a model-adequacy problem — not by an effect estimate being the wrong size.

Result: the dip at moderate thermal stress survives removal and is if anything
deeper. The rising upper limb does not survive and is a boundary artefact in both
fits. Conclusion recorded: a negative association at moderate DHW is defensible;
nothing about the high end is identified and it must not be quoted.

**2026-09-06 — Step 3 sensitivity suite specified in advance of running it.**
Five checks agreed before any were run: leave-one-region-year-out influence,
restriction to the sheltered stratum, each region fitted alone, a lagged cyclone
term, and a Tweedie refit on the density scale. Specifying the full set first
means none of them can be presented as "the one that worked".

**2026-09-06 — Cyclone integration downgraded to a lag, with the reason recorded.**
The Stage 1 plan asked for cyclone exposure integrated across survey intervals.
That is not possible from this extract: exposure is recorded only at survey points
and the intervening years are absent from the file. A lagged term (the previous
survey at the same site) is what the data support, and it is a weaker substitute.
Proper integration needs external cyclone track data, which is now on the list of
questions for AIMS.

**2026-09-06 — Habitat covariates NOT swapped despite the null result.**
`LCC_%` has more between-region-year variance than the chosen `LHC_%` and fell
harder in 2017, so swapping them was considered. Rejected: the habitat block came
out null, and changing a covariate after seeing a null result is selection on the
outcome. The null is reported as conditional on the covariates chosen in advance,
and that conditionality is stated as a limitation rather than engineered away.
