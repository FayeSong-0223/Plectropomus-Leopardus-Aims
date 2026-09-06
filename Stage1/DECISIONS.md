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
