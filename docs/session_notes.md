# Session Notes (append-only)

## 2026-01-08
- Started repo scaffolding plan for Giotto ST pipeline + Singularity container.
- Next: create skeleton directories/files; add minimal runnable “hello pipeline” script inside container.

## 2025-12-11
----------
- Environment: R 4.4.3 with user library at ~/R/4.4.3_lib; Giotto suite installed.
- Added Giotto Miniforge environment at ~/.local/share/r-miniconda/envs/giotto_env (Python 3.6) with leidenalg, python-igraph, scikit-learn, etc.
- Batch pipeline in Giotto.R processes all Xenium output-* folders, writing results to giotto_runs/<run_id>/.
- Normalization, PCA (FactoMineR backend), UMAP, nearest neighbor graph, and Leiden clustering run end-to-end per dataset.
- Outputs per run: giotto_object.rds, sessionInfo.txt, spatPlot2D.png, UMAP.png; clustering metadata (leiden_clus) verified in saved objects.
- Outstanding next steps:
  1. Review spatial and UMAP plots in giotto_runs/ for QC.
  2. Adjust filtering/normalization thresholds if low-quality cells should be removed.
  3. Add future analyses after the PCA/UMAP block in run_xenium_analysis(), writing outputs under each run directory.
- Reminders: prepend ~/R/4.4.3_lib to .libPaths(); reuse python_path = "~/.local/share/r-miniconda/envs/giotto_env/bin/python" when creating instructions; run batch via Rscript Giotto.R or individual runs via run_xenium_analysis(<h5 path>).

## 2026-01-09
- Ported legacy Xenium ingest and baseline Giotto workflow into modular helpers R/ingest_xenium.R and R/pipeline_basic.R.
- Updated scripts/run_all.R to orchestrate Xenium runs end-to-end, emit run metadata, cluster tables, Giotto object, and spatial/UMAP plots under --output_dir.
- Added deterministic plot filenames and recorded key output paths in metadata for downstream reproducibility.
- Updated Copilot instructions for `run_all.R` to treat Xenium as a first-class input format
  (distinct from Visium), using `R/Giotto.R` as the reference implementation.
- Loaded udunits/2.2.28 and hdf5/1.12.2 modules, then installed remaining Giotto dependencies (hdf5r, zigg, FactoMineR, RcppAnnoy, labeling, etc.) into ~/R/4.4.3_lib for consistent Xenium runs.
- Swapped spatial/UMAP plotting calls in R/pipeline_basic.R to Giotto::spatPlot2D/plotUMAP with save_param metadata for current Giotto API compatibility.
- Executed scripts/run_all.R on output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834, producing results/xenium_r1/ with plots, cluster table, Giotto object, and session info.
- Outstanding: install python-igraph and leidenalg inside ~/.local/share/r-miniconda/envs/giotto_env to silence Giotto Python warnings during future runs.
- Added QC metric/plot generation in R/pipeline_basic.R (histograms, scatter, metrics CSV, interpretive text) and re-ran scripts/run_all.R to populate results/xenium_r1/qc/.
- Added repo hygiene (.gitignore, removed macOS artefacts) plus initial README/docs/images scaffold and configs/README.md to mirror CAF documentation style.
- Initialized renv with project-specific package snapshot and mirrored user library to generate renv.lock; authored container/Dockerfile with helper scripts (build/run) for remote Docker builds and documented workflow.

## 2026-01-12
- Loaded site R module to regain Rscript access inside workspace shell.
- Added R/ingest_visium_hd.R to create Giotto objects from Visium/Visium HD Spaceranger outputs (supports H5 or MTX matrices, CSV or parquet tissue positions, optional images).
- Extended scripts/run_all.R auto-detection, dry-run validation, and ingest routing to call the Visium helper while keeping Xenium logic intact.
- Updated README with Visium usage instructions and noted optional arrow dependency for parquet metadata.
- Installed arrow (plus assertthat, bit, bit64) via renv using R 4.4.3 and refreshed renv.lock to capture the new dependency.
- Ran Visium HD dry-run + execute cycle on square_008um and square_016um bins; ingest completes but full Giotto pipeline is OOM-killed on login node (≈137k bins). Recommend re-running on a compute node or downsampling before normalization.
- Added --max_cells flag to downsample within Giotto; successfully executed Visium HD square_016um with --max_cells 4000, producing outputs under results/visium_hd_square_016um_max4k/.
- Re-ran Visium HD square_008um with --max_cells 6000; outputs written to results/visium_hd_square_008um_max6k/ without OOM.

## 2026-01-13
- Removed incomplete Visium HD run directories (results/visium_hd_square_008um, results/visium_hd_square_016um, results/visium_hd_square_008um_max20k) to reduce clutter.
- Replaced synthetic docs/images assets with plots exported from results/visium_hd_square_008um_max6k/ (spatial, UMAP, and QC panels) and updated README Example Outputs accordingly.
- Added h5ad ingest path with layer/raw fallbacks, refactored sparse conversion helpers, validated via synthetic AnnData dry-run, and extended README flags/usage to describe the new format.
- Introduced optional QC filtering (`--min_genes_per_cell`) before normalization, wiring CLI flag through run metadata and README documentation.
- Expanded QC thresholds to include `--min_total_expr_per_cell`, `--max_mito_pct`, and configurable `--mito_gene_prefixes`, with run metadata summarizing removals and mitochondrial gene coverage.
- Patched Giotto expression access to use the raw_exprs slot (get_expression_values removed upstream) and re-ran scripts/run_all.R on XETG00202_R1 with filters (min genes 10, min total expr 80, max mito 25), producing results/xenium_qc_filtered/ with 22,494 retained cells and populated QC summaries.
- Added per-threshold QC logging in scripts/run_all.R, exporting metadata/XETG00202_R1_qc_qc_filter_summary.csv with before/after counts for each filter and recording the path in run_parameters.json.
- Updated QC summary export filename to metadata/XETG00202_R1_qc_filter_summary.csv to avoid duplicated suffix and removed the legacy file.
- Enabled Visium ingest in scripts/run_all.R, adding ZIP-aware spatial handling in R/ingest_visium_hd.R to auto-extract archives, then executed runs on SCAF3120_23001460_Veh_IaG_4 and SCAF3121_23001459_YAP_5Wk_3 (max_cells=6000) producing results/visium_SCAF3120/ and results/visium_SCAF3121/.
- Refreshed README and QUICKSTART to highlight spatial ZIP extraction, per-filter QC summaries, and the new Visium datasets.
- Clarified setup docs to include cloning the repository before running scripts.

## 2026-04-02
- Pulled the latest GitHub commits onto local main (through f88aa5b) to sync the workspace before review.
- Reviewed the new Visium and h5ad ingest paths plus container scaffolding; confirmed Docker CLI is installed locally, but the Docker daemon was not running, so image build could not complete.
- Fixed a real container build blocker by changing container/Dockerfile to copy renv/settings.json instead of the non-existent renv/settings.dcf.
- Updated scripts/run_all.R CLI help text so --input_format advertises h5ad alongside xenium, visium, and matrix.
- Smoke-tested Rscript scripts/run_all.R --help; the entrypoint currently fails in this local environment until renv-installed packages (notably data.table) are restored.
- Aligned the repo with the multi-gene-correlations CAF pattern: disabled automatic renv activation in .Rprofile, added a .dockerignore, added container/run_apptainer.sh, and documented the canonical /data and /output bind layout.
- Updated container/build.sh so extra Docker arguments (for example --platform linux/amd64) can be forwarded from the helper script.
- Refactored the workflow controller to support `--stage validate|ingest|qc|analyze|export|all`, added `--input_object` for restartable later-stage runs, and documented the new checkpoint objects in README, Quickstart, architecture, and the decision log.
- Diagnosed Docker HTTPS failures as NIH TLS interception inside Linux containers (`NIH-DPKI-NS-SSLCA-1A` chaining to `NIH-DPKI-ROOT-1A`) rather than an R-specific bug.
- Updated container/Dockerfile to accept an optional BuildKit secret (`extra_ca`) so enterprise root certificates can be injected at build time without committing them into the repository.
- Added a configurable `RENV_RESTORE_EXCLUDE` Docker build argument (defaulting to `arrow`) so the default image can build within Docker Desktop disk limits while preserving CSV-based Visium support.
- Successfully built `giotto-st-pipeline:dev` with `./container/build.sh giotto-st-pipeline:dev --secret id=extra_ca,src=/tmp/nih-dpki-root-1a.pem`.
- Smoke-tested the built image with `docker run --rm giotto-st-pipeline:dev --help` and stage-aware error paths for `--stage analyze` and `--stage validate`; both exercised the new staged CLI and wrote run metadata under `/output/metadata/run_parameters.json` before failing on controlled missing-input cases.
- Decided to treat the first GHCR publication as a lean image that excludes `arrow` and therefore expects non-parquet Visium metadata at runtime.
- Updated public docs to describe enterprise CA handling as a maintainer rebuild concern only, using the existing BuildKit secret path in container/Dockerfile and container/README.md.
- Added OCI image labels to container/Dockerfile and removed explicit cross-repo references from user-facing docs so the repo keeps the same operational shape without naming other projects.
- Added container/publish.sh plus GHCR publication docs describing the lean-image tag convention: one immutable `sha-<git-sha>` tag and floating `lean` plus `latest` aliases.
- Hardened container/publish.sh to resolve the local source image by image ID instead of relying on direct tag inspection, which was inconsistent on the current Docker Desktop daemon.

## 2026-04-08
- Confirmed `ghcr.io/nidap-community/giotto-st-pipeline:latest` is published and runnable via `docker pull` plus `docker run --rm ... --help`.
- Updated README and QUICKSTART to make the published GHCR image the default pull/run path, with Docker and Singularity examples that use the `/data` and `/output` contract.
- Updated container/README.md so normal users pull the published image, while source rebuild and GHCR publication guidance remain maintainer-facing.
- Recorded the GHCR publication milestone in CHANGELOG.md.
- Tested the published image against the small h5ad fixture and found a real container defect: `anndata` is missing from the image's Python environment, so h5ad ingest fails even though CLI help succeeds.
- Patched the Dockerfile to create an in-image Python virtual environment at `/opt/giotto-python` with `anndata` and `scipy`, and set `RETICULATE_PYTHON` to that interpreter for local rebuild testing.
- Found a second rebuild issue after invalidating the cached apt layer: the base image was still using plain HTTP Ubuntu mirrors, which fail signature validation on this network path. Patched the Dockerfile to rewrite Ubuntu mirrors to HTTPS before `apt-get update`.
- Found a third rebuild issue: the optional enterprise CA was being installed after `apt-get`, so it could not help HTTPS mirror access. Reordered the Dockerfile so CA injection happens before package installation.
- Switched the default PCA backend from `FactoMineR` to `irlba` in the active pipeline code so the lean analysis path no longer depends on the heavier `FactoMineR` stack.
- Added a manual GitHub Actions workflow at `.github/workflows/publish-ghcr.yml` so the lean image can be rebuilt and republished on GitHub-hosted runners when local Docker rebuilds are blocked by enterprise network issues.
- Extended the Dockerfile's baked Python environment to include the Giotto-related runtime modules needed for `h5ad` ingest and clustering (`anndata`, `scipy`, `python-igraph`, `leidenalg`, `networkx`, `scikit-learn`, `python-louvain`) and added `cmake` for more reliable package builds in GitHub-hosted rebuilds.
- The first GitHub-hosted publish attempt failed with `permission_denied: write_package` while using `GITHUB_TOKEN`, so the workflow was updated to prefer an Actions secret named `GHCR_TOKEN` for GHCR pushes and fall back to `GITHUB_TOKEN` only when no PAT is configured.
- Added `libgmp3-dev` and `pandoc` to the Dockerfile after reviewing the build-failure analysis, since the lockfile includes `gmp` and `rmarkdown` and those system dependencies may be needed during clean GitHub-hosted rebuilds.
