## 2026-01-12
- Added Visium/Visium HD ingest path in scripts/run_all.R with supporting helper for Spaceranger outs.
- Introduced optional --max_cells downsampling and validated Visium HD run on square_016um (4k spots) dataset.

## 2026-01-13
- Added per-threshold QC summary export (`<project_id>_filter_summary.csv`) and logged filter provenance in run metadata.
- Extended Visium ingest to auto-extract Spaceranger `*spatial*.zip` assets; executed SCAF3120 and SCAF3121 datasets from GLOBUS import with downsampling support.
