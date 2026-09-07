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

**2026-09-06 — Overreach recorded: a lag term was fitted that the plan had cut.**
Section 13 of the project plan cuts "lag structures and distributed-lag models — 13
region-years cannot identify them." A lagged cyclone term was fitted in Step 3 anyway,
justified at the time as a substitute for the interval integration that turned out to
be impossible. That justification does not hold: the impossibility of integration does
not make a lag identifiable, and the plan had already ruled it out on those grounds.
The term was not significant (p = 0.149) so nothing rests on it, but it should have
been flagged as crossing a stated limit rather than presented as a sensible substitute.

**2026-09-06 — Four outstanding plan deliverables completed in Step 4.**
An audit against the project plan found four specified items had not been produced:
the effect magnitudes Q1 actually asks for, the Moran's I spatial residual check
(section 10, check 4), the variable dictionary (Table 1), and the shared components
of the variance decomposition in the figure rather than only the log. All four are
now in coral_trout_stage2_step4.R. None required new data.

**2026-09-06 — Spatial structure found in Whitsunday residuals; NOT acted on here.**
The Moran's I check the plan asked for found residual spatial autocorrelation in
Whitsunday (I = 0.146, p = 0.007) but not in Palm (I = 0.011, p = 0.878). The site
random effect has therefore not absorbed the spatial signal in one of the two regions,
and a spatial term would be justified there. This is left as a documented finding
rather than fixed, because adding a spatial field is the first item section 13 cuts
and doing it now would be a substantial scope change made at the end of the project.
It goes to the top of any version two.

**2026-09-06 — Full review found five things, three of them real problems.**

*Corrected.* The Moran's I result was seed-dependent. Randomised quantile residuals draw
from a uniform, every script fixed `set.seed(1)`, and no second randomisation was ever
checked. Across 20 seeds Whitsunday's I ranges 0.079–0.167 and the permutation p falls
below 0.05 in 11 of 20; seed 1 sat near the top. The claim has been downgraded from
"spatial structure detected" to "consistently positive, not established". Palm is near
zero and never significant under any seed, so the *contrast* between regions holds.

*Corrected.* `s(depth)` had concurvity 1.000 with the site random effect, because depth
is constant within every site in this subset. Depth is now a linear term: identical
estimates, AIC 2353.0 vs 2353.1, no aliasing — the smooth had already collapsed to a
straight line (edf 1.005). Omitting depth raises the Whitsunday estimate to 3.62, so it
is doing real adjustment work and is kept.

*Corrected.* The report quoted Whitsunday rugosity moving 3.63 to 3.77 between 2016 and
2017. The correct values are 3.58 and 3.72. The point being made — that structural
complexity did not fall while coral cover halved — is unaffected.

*Newly reported, not previously computed.* Concurvity was never checked. Worst-case
values are 0.951 for the thermal smooth and 0.926 for turbidity against the
region-specific year trends, which quantifies a confounding that had only been described
qualitatively.

*Newly reported, not previously computed.* Basis dimension was never checked. k.check
flags the year smooths and kd490. For the year smooths this cannot be fixed by raising k
(Palm has six distinct years) and Model B is the specification that answers it; for
kd490 it is a genuine limitation.

Verified and unchanged: every headline figure was recomputed independently from the raw
CSVs by a separate code path, and all 30-plus matched. A clean end-to-end run of all six
scripts reproduces all 20 tables and logs byte for byte, with no warnings.

---

## Second audit, 2026-09-06

A read-only audit of the whole project against both planning documents. Seven claims
were tested independently; five held, one held only partly, and one was wrong. What
follows is what changed as a result. The standing rule still applies: every change
below was prompted by a defect in method or a mismatch with the data, none by an
effect estimate.

**Withdrawn — the variance decomposition answering Q2.** The reported figure was
"management explains 48.1% of the variation". Three faults, any one of which is
disqualifying:

1. It was not invariant to how the factor contrasts were coded. Refitting the
   identical model under sum-to-zero rather than treatment contrasts moved management
   from 0.481 to 0.391 and space/time from 0.121 to 0.345, with identical fitted
   values and identical log-likelihood. A number that moves when nothing about the fit
   moves describes the parameterisation, not the data.
2. It partitioned the systematic part only. The linear predictor has variance 0.533;
   the negative binomial observation process contributes about 0.455 more on the same
   scale. Roughly half the variation was outside the partition and the shares were
   silently conditional on that.
3. The site block used the empirical variance of the shrunken random-effect
   predictions, 0.048, against an estimated variance component of 0.094. Site was
   understated by about half and every other block inflated against it.

The report had carried a caveat on this number, but it gave a different and less
serious reason than any of the above. `table5_variance_decomposition.csv` and the old
figure 7 are gone.

**Replaced by a drop-one-block comparison** — and that replacement was itself
superseded within a day. See the third audit below. The version described here
compared deviance explained across models that each re-estimated theta, and measured
management by removing it from a model with no random effect. Both are wrong, for
reasons set out in the next section. Nothing from this paragraph should be quoted.

**Environmental intervals were too narrow.** Every interval in the project used
`vcov(m)`, which treats the smoothing parameters as known. With
`vcov(m, unconditional = TRUE)`, rugosity (1.27, 0.98–1.64), turbidity (0.75,
0.55–1.01) and cyclone exposure (0.85, 0.67–1.06) all cease to exclude 1. Only the
thermal term survives, at 0.68 (0.48–0.96). The unconditional matrix is now used
everywhere, including for protection, where it widens the Whitsunday interval from
2.28–4.52 to 2.25–4.58 and changes nothing.

**"Survives every check" was not true of Palm.** The sheltered-only restriction gives
Palm NTR 1987 as 0.47 (0.24–0.93) — significantly *below* 1, where every other
specification gives 1.1 to 1.4. Section 10 of the plan set the rule in advance: "If it
agrees with the adjusted estimate, say so. If it does not, distrust the adjusted one."
That rule was not applied. The report noted the anomaly in its sensitivity section
while its summary still claimed the contrast survived everything, and the README's
"excludes 1?" column recorded "no" for a row that does exclude 1. Corrected in both.
The Whitsunday result is unaffected: it survives the restriction at 2.58 (1.56–4.25).

**Two published values for one model.** Step 2's protection ladder labelled its first
rung "baseline (Step 1)" but included a depth term that Step 1 did not have. Table 3
therefore gave Whitsunday NTR 1987 as 3.59 and table 4 gave 3.16 for what was named as
the same model. The ladder now starts with the Step 1 model exactly as Step 1 fits it,
and depth enters as its own rung — which also makes visible that depth moves the
estimate more than any other single covariate.

**Keppel 2021 miscounted.** The data question sent to researchers said legal density
exceeds total density in 18 rows. It is 17, at ratios of 2.1 to 6.8, all in Keppel
2021. The 18th row is Whitsunday 2018 site HY3, where legal exceeds total by 0.0018 —
a rounding artefact in a different region-year, and not part of the same phenomenon.
The comparison range is 1.10–2.17 kg, not 1.11–2.17.

**Stage 1 note rugosity.** The note still quoted Whitsunday rugosity moving 3.63 to
3.77 between 2016 and 2017. The correct values, 3.58 and 3.72, had been fixed in the
report months earlier and never propagated. Fixed.

**The environment is now recorded.** Every log carries R version, mgcv version and
platform. "No packages are required" was true and misleading: mgcv ships with R, but
its version tracks the R version, and REML fitting and the `nb()` family have both
changed across releases.

### Newly run, never previously checked

**Excess zeros.** The plan flagged Palm's concentrated zeros in section 8 and no
distributional check ever followed. Tested by parametric bootstrap from the fitted
model: 52 zeros observed against a simulated median of 35 (95% interval 25–46,
p = 0.004). Palm alone p = 0.009, Whitsunday p = 0.062. The excess is real and is a
stated limitation. It does not drive the result — a zero excess inflates apparent
overdispersion, which widens intervals rather than narrowing them, and trimming the
one site that is mostly zeros leaves Whitsunday at 3.27 (2.31–4.63) against 3.21.

**Temporal autocorrelation within sites, and it is not clean.** The model gives each
site a random intercept and nothing had tested whether that is enough. Pairs are
binned by the true gap in years, because the survey years are unequally spaced and
pairing consecutive surveys would call a one-year gap and a three-year gap the same
thing. Ten of eleven gaps sit at or above the correlation a site random intercept
induces mechanically (about −0.18 here). The three-year gap does not: r = −0.376,
significant under all 20 randomisations, carried by the 2009-to-2012 interval in both
regions. A reversal of that size means something moved sites differentially between
those surveys and the model has no term for it. This is the cyclone limitation made
concrete: exposure is recorded at survey points, so a disturbance falling inside a
three-year gap is invisible to the covariate. Reported as a limitation and added to
the questions for anyone holding track data. Naming a particular storm would be
inference this dataset cannot support.

**Multiplicity.** Model A reports 17 tests and no p-value quoted anywhere in this
project had been adjusted. One correction across all 17 would be the wrong instrument,
since those tests answer four different questions. Holm within each of the plan's four
hypothesis families instead. Q3 — the region-by-protection interaction — clears it by
an order of magnitude (adjusted p = 0.0014 and 0.0020). Nothing in the habitat family
survives. In the environment family the thermal term survives both this and the
unconditional interval; cyclone exposure survives neither; turbidity survives Holm but
not the interval, so it is reported as unresolved. Where two corrections disagree, the
weaker reading is the one quoted.

### Deliberately not changed

*The figure and table count.* Section 11 of the plan specifies six figures and two
tables and says to cut any figure that does not answer one of the four questions. The
project has sixteen figures and fourteen tables. Every addition since has been a
diagnostic or a correction rather than a new result, and removing the evidence for a
correction to make a count would be the wrong trade. Recorded as a deviation rather
than defended as compliance.

*The protection estimate itself.* Nothing about the Whitsunday result changed under
any of the above: not the unconditional intervals, not the multiplicity correction,
not the zero-count trimming, not the corrected baseline ladder. It was not adjusted to
keep it that way; it simply did not move.

### What the audit got wrong

One claim tested did not hold. Outputs were said not to match the current scripts. A
clean run from an empty `outputs/` reproduces every log, table and figure byte for
byte, so the scripts and their outputs are consistent. What was true is narrower and
was found separately: the copies committed to the repository for steps 2 and 3 predate
the depth respecification, so the repository — not the pipeline — was carrying stale
logs. Fixed by committing a complete regenerated set.

---

## Third audit, 2026-09-07

The second audit's replacement for Q2 was itself wrong, and two of its new diagnostics
were built on a simulation that did not respect the model's own structure. This round
fixes those. It also removes dead code that would have silently rebuilt a withdrawn
result. Where a claim from the previous round is overturned, it is marked superseded
above rather than deleted.

**Dead code that made the previous reproducibility claim false.** Step 4 still
contained the withdrawn variance decomposition and still wrote
`fig7_variance_decomposition.png`. The figure had been deleted by hand from the output
directory, so the directory looked correct while the pipeline that produced it did not:
any clean run would have recreated the file. The previous statement that two runs were
byte-for-byte identical was therefore not a statement about the analysis as it now
stands. The section is removed, and the reproducibility claim is only made after two
clean runs of the current code.

### Q2, rebuilt for the second time

**Fault 1: deviance explained was not comparable across the models being compared.**
`nb()` estimates theta by REML for each model separately. A reduced model can absorb
the structure it lost into a smaller theta, so its deviance is measured against a
different variance function and the difference between two deviances is not a
difference in fit. Every model in the drop-one-block comparison is now fitted with
`negbin(theta)` at the full model's estimate, so all deviances share one scale.

**Fault 2: the 14.09 percentage-point "management" figure was confounded.** It came
from removing management from a model with no site random effect. With no random effect
present, the management terms absorb every between-site difference that happens to
align with zoning — habitat, history, location, anything fixed about a site. It
measures confounding as much as management. Reporting it as management's contribution
would have repeated, in a new form, the error the second audit was meant to fix. It is
withdrawn.

**What is reported now.** Three things, in order of how much weight they carry.

*First, a statement of what is not identifiable.* Protection changes within 0 of 71
sites. Its entire contribution is between-site, and a site random intercept is also
entirely between-site; the two are estimated from the same degrees of freedom. No
partition of variance or of deviance can say how much of a between-site difference
belongs to zoning rather than to whatever else distinguishes those sites. This is a
limit of the design and cannot be fixed by a better statistic. The project does not
report a management percentage, and the honest answer to "how much does management
explain" is that this design cannot say.

*Second, drop-one-block at fixed theta*, with a site-level cluster bootstrap for
intervals. Descriptive and non-additive: the drops do not sum to anything and are not
unique variance shares. The management row is expected to sit near zero for the reason
just given, and is not evidence of no effect.

*Third, out-of-sample prediction to unseen sites*, which is the one comparison the
overlap does not spoil. Ten folds split by site, so a site is never in both training and
test; each model scored by negative binomial log predictive density on held-out sites
with the random effect set to zero, because an unseen site has no intercept to estimate.
This is a predictive statement only. It does not partition variance, it does not license
a causal reading, and it cannot be used to argue that management contributes more
variance than site identity. The site random effect is deliberately excluded from the
blocks compared, because it is already absent from every prediction and its row would
invite exactly that unsupported comparison.

### Diagnostics, rebuilt

**The parametric bootstrap did not preserve site clustering.** The previous version
simulated from `fitted(mA)`, which bakes the shrunken site estimates into every
replicate as if they were known constants. The empirical standard deviation of those
fitted effects is 0.219 against an estimated component of 0.307, so the simulated data
carried about 30% less between-site variation than the model implies, and the null
distributions were too narrow. Every replicate now draws a fresh set of site effects
from N(0, sigma_site^2), one per site, so clustering is generated rather than
inherited. The correction is not cosmetic: it widened the null for the residual
correlation statistic from roughly (−0.21, −0.15) to (−0.32, −0.20).

**The zero check was a conditional simulation, not a parametric bootstrap.** It compared
the observed zero count against draws at fixed mu, which ignores that mu and theta were
estimated from the same data. Each replicate now simulates, refits, and recomputes the
statistic from the refit exactly as it is computed from the real fit.

**The temporal check searched eleven gaps without adjustment and treated jitter seeds as
evidence.** Correlations were tested against zero, when a site random intercept induces
a negative correlation at every gap by construction; the strongest of eleven gaps was
selected after seeing all of them; and twenty randomised-quantile seeds were reported as
though they were twenty independent tests, when they are twenty draws from one dataset.
The statistic is now the most negative correlation across all gaps — so the search is
part of the statistic — compared against its own simulated null, which carries both the
induced correlation and the search. One jitter draw per replicate, on both sides.

**No mechanism is proposed for what the temporal statistic shows.** The previous version
described a disturbance falling inside a survey gap and connected it to the cyclone
covariate's limitations. That is an ecological explanation this analysis does not test:
there is no covariate here that would distinguish it from any other source of the same
pattern. The finding is reported as what it is, a statement about whether the site
random intercept is sufficient, and the survey pairs are listed descriptively.

**Simulated p-values.** Computed as (extreme + 1) / (B + 1) and never reported as zero;
the smallest attainable value is stated alongside. Where a p-value falls between 0.01
and 0.10 the result is labelled uncertain rather than pushed across a threshold. More
replicates are not added to resolve such a case unless the answer would change a
conclusion.

**Planned sensitivity now run.** A site-level random slope on year, so that sites are
allowed their own trajectory rather than only their own level, alongside the
zero-trimmed refit.

### Other corrections

**Basis dimension is tested rather than asserted.** The previous version stated that the
year-smooth flag could not be fixed by raising k and that the kd490 flag was a real
limitation. Neither was checked. Both are now checked by refitting at the largest basis
the data allow — k up to 20 for kd490, k up to the number of distinct survey years for
the year smooths — and reporting whether the flag clears and whether the protection
estimate moves.

**Observed against fitted for the final model.** Step 1 checked this for the baseline;
the final model had never been checked the same way. Binned means with intervals.

**Partial-effect bands were conditional while the text claimed unconditional
everywhere.** Figures 6 and 9 drew `plot.gam` bands from the conditional covariance
matrix. Both now pass `unconditional = TRUE`, matching the intervals in the tables.

**The Keppel anomaly count is now reported in full.** There are 18 rows in the extract
where legal density exceeds total density: 17 in Keppel 2021 at factors of 2.1 to 6.8,
and one in Whitsunday 2018 where the excess is 0.0018 — a rounding artefact, three
orders of magnitude smaller and a different thing entirely. The question sent to AIMS
names 17, and the code now reports both figures so that anyone repeating the check
finds 18 and is not left wondering which is right.

**Language.** Whitsunday is a robust adjusted association, not a measured effect of
protection: adjusting for wave exposure and restricting to the sheltered stratum rule
out exposure as the explanation, but reserve placement was not random and no
observational check here establishes that whatever else distinguishes the zoned sites is
innocuous. Palm is unresolved, not null. Q2 yields no identifiable management
percentage.
