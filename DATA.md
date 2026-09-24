# Data provenance and licence

## Source

The data come from the Australian Institute of Marine Science (AIMS) dataset
**Spatio-temporal dynamics of coral reef fish assemblages on inshore reefs of the Great Barrier Reef**. Its metadata record is:

<https://apps.aims.gov.au/metadata/view/814a0be3-ed85-4a43-87b7-59ece6eb6a05>

Both files below were downloaded from that record's *Data Downloads* section on
2 September 2026. This analysis uses the Palm and Whitsunday surveys from 2007 to 2018.

| File | Contents |
|---|---|
| `Selected fish benthic physical sitelevel 2021.csv` | Site-level fish, benthic and physical variables. 637 rows, 304 columns. |
| `Inshore fish site coordinates.csv` | Site coordinates in degrees-decimal-minutes. 112 rows. |

**Neither file is included in this repository.** To reproduce the analysis, place both
in `Stage1/data/` and run the scripts from the `Stage1/` directory, starting with
`coral_trout_stage1.R`.

**Citation.** As the record asks:

> Australian Institute of Marine Science (AIMS). (2022). Spatio-temporal dynamics of coral reef fish assemblages on inshore reefs of the Great Barrier Reef. https://apps.aims.gov.au/metadata/view/814a0be3-ed85-4a43-87b7-59ece6eb6a05, accessed 02-Sep-2026.

The record also lists two papers that describe the survey programme and analyse
these data:

- Ceccarelli, D. M., Evans, R. D., Logan, M., Jones, G. P., Puotinen, M., Petus, C.,
  Russ, G. R., Srinivasan, M., & Williamson, D. H. (2023). Physical, biological and
  anthropogenic drivers of spatial patterns of coral reef fish assemblages at regional
  and local scales. *Science of The Total Environment*, 904, 166695.
  https://doi.org/10.1016/j.scitotenv.2023.166695
- Ceccarelli, D. M., Logan, M., Evans, R. D., Jones, G. P., Puotinen, M., Petus, C.,
  Russ, G. R., Sinclair-Taylor, T., Srinivasan, M., & Williamson, D. H. (2024).
  Regional-scale disturbances drive long-term decline of inshore coral reef fish
  assemblages in the Great Barrier Reef Marine Park. *Global Change Biology*, 30,
  e17506. https://doi.org/10.1111/gcb.17506

**Data provenance.** The source extract and associated site-coordinate file used in this project were downloaded directly from publicly available Australian Institute of Marine Science (AIMS) data resources in 2026. No restricted, confidential or privately supplied AIMS data are used in this repository. The source data remain subject to the licence and attribution requirements specified by AIMS for the original dataset.

## Licence

The source data are published by AIMS and are not redistributed here.

AIMS sets the licence for each dataset in its metadata record. The record for this
dataset lists it under *Constraints* as "Attribution 3 Australia", that is the
**Creative Commons Attribution 3.0 Australia (CC BY 3.0 AU)** licence.

## Attribution

Material in `Stage1/outputs/` — the figures, summary tables and logs — is derived
from AIMS data. AIMS requires the following attribution for modified or derived
material:

> Based on Australian Institute of Marine Science data.

Please carry this attribution, and cite the dataset as above, if you reuse the figures
or tables.

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
The derived figures and tables remain subject to the dataset's CC BY 3.0 AU licence
and the attribution described above.

CC BY carries no share-alike requirement, so original work here may be licensed
freely. It does not, however, permit granting others broader rights to the underlying
AIMS material than are held in the first place.

## Contact

Questions about the source data: adc@aims.gov.au
