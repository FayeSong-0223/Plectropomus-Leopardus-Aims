# Data provenance and licence

## Source

This analysis uses an extract from the Australian Institute of Marine Science (AIMS)
inshore reef monitoring programme, covering 2007–2018:

| File | Contents |
|---|---|
| `Selected fish benthic physical sitelevel 2021.csv` | Site-level fish, benthic and physical variables. 637 rows, 304 columns. |
| `Inshore fish site coordinates.csv` | Site coordinates in degrees-decimal-minutes. 112 rows. |

**Neither file is included in this repository.** To reproduce the analysis, place both
in `Stage1/data/` and run the scripts from the `Stage1/` directory, starting with
`coral_trout_stage1.R`.

**Data provenance.** The source extract and associated site-coordinate file used in this project were downloaded directly from publicly available Australian Institute of Marine Science (AIMS) data resources in 2026. No restricted, confidential or privately supplied AIMS data are used in this repository. The source data remain subject to the licence and attribution requirements specified by AIMS for the original dataset.

## Licence

The source data are © Australian Institute of Marine Science and are not redistributed
here.

AIMS licenses each dataset individually in its metadata record, and the licences are
not uniform. General AIMS website material is CC-BY 4.0 Australia, but the reef
monitoring fish series from which this extract descends is published under a
**Creative Commons Attribution–NonCommercial (CC-BY-NC) 4.0** licence, and some
inshore Marine Monitoring Program records carry further access and use constraints.
AIMS also notes that additional terms beyond the Creative Commons licence may apply.

## Attribution

Material in `Stage1/outputs/` — the figures, summary tables and logs — is derived
from AIMS data. AIMS requires the following attribution for modified or derived
material:

> Based on Australian Institute of Marine Science data.

Please carry this attribution if you reuse the figures or tables.

## What is published, and what is not

**Published.** The six analysis scripts, sixteen figures, the summary tables in
`Stage1/outputs/tables.md`, the Stage 1 audit log and the five Stage 2 console logs.
These are results rather than data: every table is an aggregate of about twenty rows
or fewer, and the logs contain only summary statistics and model output.

**Not published.** The two source CSVs, and `Stage1/outputs/analysis_dataset.csv`.

That last file is 467 rows by 20 columns — site identifiers, coordinates to five
decimal places, every covariate, and the recovered integer counts. Every column is
either copied from the AIMS extract or a deterministic transform of it, so it is a
filtered copy of the source rather than a summary, and publishing it would amount to
redistribution. Withholding it costs little: the script regenerates it in one command
from the source files.

## Licensing of this repository

The R code and the written text are released under the MIT Licence (see `LICENSE`).
The derived figures and tables remain subject to the AIMS terms described above.

CC-BY-NC carries no share-alike requirement, so original work here may be licensed
freely. It does not, however, permit granting others broader rights to the underlying
AIMS material than are held in the first place.

## Contact

Questions about the source data: adc@aims.gov.au
