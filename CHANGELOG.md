## 2026-06-13
- Reframed user instructions around local command-line R and R Markdown first, followed by Singularity/Apptainer and Docker container paths.
- Added reusable YAML config templates for Xenium, Visium/Visium HD, and h5ad runs, plus an example R Markdown workflow.
- Added automatic publication provenance metadata at `metadata/provenance.json` and embedded the same provenance block in `metadata/run_parameters.json`.
- Updated Apptainer/Singularity helper behavior to pass runtime, SIF path, SIF SHA-256, and optional source image digest into run metadata.

## 2026-06-12
- Split documentation into a user-facing root README, quickstart, container guide, and developer README with updated example output figures.
- Added CLI controls for plot sizing and core analysis settings: `--spatial_point_size`, `--umap_point_size`, `--spatial_legend_text`, `--spatial_axis_title`, `--pca_dims`, `--neighbor_k`, and `--cluster_resolution`.
- Improved Visium ingest by honoring `in_tissue` spot filtering, dropping all-zero genes after barcode matching, and preserving restartable stage outputs.
- Fixed Giotto UMAP execution by passing an explicit thread count into `runUMAP()`.
- Made container builds more reproducible by pinning the Rocker base image digest, restoring the full `renv.lock` by default, pinning the `renv` bootstrap version, and adding pinned Python requirement files.
- Added OCI provenance labels for source commit and `renv.lock` SHA-256, plus GHCR publishing guards that refuse dirty local worktrees and verify local image metadata before pushing.
- Added `container/export_tls_chain_ca.sh` to generate a temporary local CA bundle for Docker builds behind TLS inspection proxies, including Bioconductor redirect targets.
- Updated GHCR tagging guidance to use immutable `sha-<full-git-sha>` tags and optional registry digest pins for reproducible runs.
- Validated the real Visium smoke workflow and export-only rerun against `/private/tmp/giotto-readme-v1-example`; local Docker Bioconductor SSL access now passes with the generated extra CA PEM, and a lean `arrow`-excluded Docker smoke image builds and loads `Giotto` successfully.

## 2026-01-12
- Added Visium/Visium HD ingest path in scripts/run_all.R with supporting helper for Spaceranger outs.
- Introduced optional --max_cells downsampling and validated Visium HD run on square_016um (4k spots) dataset.

## 2026-01-13
- Added per-threshold QC summary export (`<project_id>_filter_summary.csv`) and logged filter provenance in run metadata.
- Extended Visium ingest to auto-extract Spaceranger `*spatial*.zip` assets; executed SCAF3120 and SCAF3121 datasets from GLOBUS import with downsampling support.

## 2026-04-08
- Published the lean GHCR image at `ghcr.io/nidap-community/giotto-st-pipeline:latest` with OCI metadata and traceable tag aliases.
- Updated README and QUICKSTART to use pull/run-first Docker and Apptainer/Singularity workflows, while keeping source rebuild guidance in `container/README.md`.
