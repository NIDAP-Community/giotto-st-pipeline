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
