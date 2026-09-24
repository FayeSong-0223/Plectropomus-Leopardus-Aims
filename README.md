# Coral trout study — *Plectropomus leopardus*, inshore Great Barrier Reef

Ziqi (Faye) Song

Analysis of AIMS inshore reef monitoring data for *Plectropomus leopardus* on the Palm
and Whitsunday inshore reefs, 2007–2018.

## Main result

On Whitsunday inshore reefs, sites inside no-take zones hold roughly three times the
density of *P. leopardus* found on fished sites. In the Palm region no difference is
detectable, and the specifications disagree. This is a robust association, not a
measured effect of protection: reserve placement was not random.

The full results, including what the analysis could not establish, are in the
**[project write-up](Stage1/README.md)**.

## Stages

Both stages are complete. All scripts, outputs and the decisions log are in
[`Stage1/`](Stage1/); the folder keeps its original name.

- **Stage 1 — data preparation and design audit** (`coral_trout_stage1.R`).
  Establishes what the data can and cannot support, audited before any model is
  fitted, and writes the analysis dataset that Stage 2 uses. Base R. Figures 1–3.
- **Stage 2 — hierarchical count model** (`coral_trout_stage2_step1.R` to
  `coral_trout_stage2_step5.R`). A negative binomial mixed model with a site random
  intercept, fitted by REML with `mgcv`, followed by habitat and environment terms,
  influence and sensitivity checks, effect magnitudes and spatial checks, and a final
  round of diagnostics. Figures 4–16.

Every change to the analysis plan, and what prompted it, is logged in
[`Stage1/DECISIONS.md`](Stage1/DECISIONS.md).

## Reproducing

The analysis needs only base R and `mgcv`, which ships with R. Place the two source
CSVs in `Stage1/data/` (see [DATA.md](DATA.md)), then run from `Stage1/`:

```bash
Rscript coral_trout_stage1.R        # run first: writes the analysis dataset
Rscript coral_trout_stage2_step1.R
Rscript coral_trout_stage2_step2.R
Rscript coral_trout_stage2_step3.R
Rscript coral_trout_stage2_step4.R
Rscript coral_trout_stage2_step5.R
```

Steps 2 and 5 are slow, roughly 50 minutes together on two cores; the rest take
seconds.

## Data

The data come from the AIMS dataset *Spatio-temporal dynamics of coral reef fish assemblages on inshore reefs of the Great Barrier Reef*
([metadata record](https://apps.aims.gov.au/metadata/view/814a0be3-ed85-4a43-87b7-59ece6eb6a05)).
The source files are **not included in this repository**. See **[DATA.md](DATA.md)**
for provenance, licence, attribution, and how to place the files to reproduce the
analysis.

Derived figures and tables are based on Australian Institute of Marine Science data.

## Licence

Code and text: MIT, see [LICENSE](LICENSE). Derived outputs remain subject to the
dataset's licence (CC BY 3.0 AU) and the attribution set out in [DATA.md](DATA.md).
