# *Plectropomus leopardus* on inshore reefs of the Palm and Whitsunday regions

**Stage 1 — data preparation and design audit.** A preparatory analysis of the AIMS
inshore reef monitoring extract, 2007–2018.

Ziqi (Faye) Song · September 2026

---

## What this is

A complete, reproducible audit of a monitoring dataset, carried out *before* any model
is fitted. It establishes what the data can and cannot support, and surfaces the design
problems that would otherwise be discovered halfway through an analysis.

Nothing here is causal. Patterns are described, not explained.

## Scope

- **Species:** *Plectropomus leopardus* only. Not pooled *Plectropomus*.
- **Regions:** Palm and Whitsunday. Keppel and Magnetic are excluded — both are almost
  entirely *P. maculatus*, so pooling would conflate species with region.
- **Response:** density (rescaled counts). Species-level biomass and size structure do
  not exist in this extract; they are available only pooled across all three congeners.
- **Records:** 467 site-years, 71 sites, 13 region-years. The panel is complete.

## Main findings from the audit

1. **Densities are rescaled counts, and the scaling constant is not uniform.** It is
   0.6666 in 2007–2017 and 0.6660 in 2018. Deriving the step per year recovers integers
   exactly (max deviation 0.00000); assuming a single step corrupts one year silently.
   Survey area is constant, so a later count model needs no offset.

2. **Protection is confounded with wave exposure, unevenly by region.** In Whitsunday,
   14 of 20 fished sites are exposed while 10 of 12 long-standing reserve sites are
   sheltered. In Palm the levels are far more comparable. The region with the large
   apparent protection effect is the one where the comparison is least trustworthy.

3. **Thermal covariates are near-perfectly nested within region-years.** `maxDHW` and
   `SSTmean` carry 99% and 98% of their variance between region-years, so they rest on
   roughly a dozen contrasts rather than 467 observations. `Cyclone` does not (18%) and
   is much better identified.

4. **Coral cover and structural complexity decoupled in 2017.** Whitsunday live hard
   coral fell from 42% to 19% while mean rugosity did not decline. Covariate choice
   between the two is therefore consequential, not cosmetic.

## Files

| File | Contents |
|---|---|
| `coral_trout_stage1.R` | The whole analysis. Base R only — no package dependencies. |
| `data/` | The two source CSVs. **Not included in this repository** — see [DATA.md](../DATA.md). |
| `outputs/fig1_design_map.png` | Sites by protection status and exposure. |
| `outputs/fig2_density_protection.png` | Density over time by region and protection. |
| `outputs/fig3_habitat_disturbance.png` | Habitat condition and disturbance exposure. |
| `outputs/table1_design_summary.csv` | Sites and observations by region, protection, exposure. |
| `outputs/table2_covariate_variance_structure.csv` | Between-region-year variance share per covariate. |
| `outputs/analysis_dataset.csv` | Analysis-ready dataset with recovered counts. **Not included** — the script regenerates it; see [DATA.md](../DATA.md). |
| `outputs/audit_log.txt` | Full console log — every number quoted above. |

## Reproducing

Place the two source CSVs in `data/`, then:

```bash
Rscript coral_trout_stage1.R
```

Run from this directory (`Stage1/`). No packages are required, so this works on a fresh R
installation. Everything in `outputs/` is regenerated.

Base R was used deliberately: it makes the script portable and immune to package
version drift. Stage 2 will need `glmmTMB` or `mgcv`, and `sf` for a coastline map.

## Data preparation issues encountered

Recorded because each is a silent failure if missed:

- The coordinate file carries a UTF-8 byte-order mark on its header, but the degree
  symbols in the coordinate strings are not valid UTF-8. Reading the file as UTF-8
  fails on the body; reading it plainly leaves the first column named `<BOM>REGION`
  so that `coord$REGION` returns `NULL` without error.
- The coordinate file contains `WHITUSNDAY` alongside `WHITSUNDAY`.
- The `NTR` field contains a double-space variant of `NTR 1987`, creating a spurious
  fourth protection level.
- The density scaling constant changes between 2017 and 2018.

## A note on the wider file (out of scope)

For Keppel 2021, `Plectropomus legal density` exceeds `Plectropomus total density` in
18 rows, and implied mean weight of legal-sized fish is 0.30 kg against 1.11–2.17 kg in
every other region-year. Total density for the same rows recovers cleanly, so the survey
area is unchanged, and the discrepancy is not a constant multiple — a units error does
not fit. Reported as a question for AIMS, not as a correction.

## Next step

A hierarchical count model with site-level random effects, protection and exposure
entered together, and honest treatment of the fact that thermal covariates rest on
about a dozen contrasts. The audit above determines what that model may claim.

---

Based on Australian Institute of Marine Science data. See [DATA.md](../DATA.md) for
provenance, licence and attribution.
