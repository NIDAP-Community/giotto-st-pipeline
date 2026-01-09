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
