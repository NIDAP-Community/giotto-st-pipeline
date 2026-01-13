# Quickstart

1. Load R and restore dependencies (once per environment):

	```bash
	module load R/4.4.3
	cd /path/to/giotto-st-pipeline
	Rscript -e 'renv::restore(prompt = FALSE)'
	```

2. Run a Xenium dataset (
	`cell_feature_matrix.h5` + `cells.csv[.gz]`) and capture outputs in a mounted directory:

	```bash
	Rscript scripts/run_all.R \
	  --input_format xenium \
	  --input_dir /data/xenium/output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834 \
	  --output_dir /data/results/xenium_r1 \
	  --project_id XETG00202_R1
	```

3. Run a Visium / Visium HD Spaceranger outs directory (the pipeline auto-unzips `*spatial*.zip` bundles if present). Downsample large runs with `--max_cells`:

	```bash
	Rscript scripts/run_all.R \
	  --input_format visium \
	  --input_dir /data/visium/SCAF3120_23001460_Veh_IaG_4/outs \
	  --output_dir /data/results/visium_SCAF3120 \
	  --project_id SCAF3120_23001460_Veh_IaG_4 \
	  --max_cells 6000
	```

4. Run an AnnData file (requires `anndata` + `scipy` in the Giotto Python environment):

	```bash
	Rscript scripts/run_all.R \
	  --input_format h5ad \
	  --input_path /data/visium/sample.h5ad \
	  --output_dir /data/results/sample_h5ad \
	  --project_id sample_h5ad \
	  --python_path ~/.local/share/r-miniconda/envs/giotto_env/bin/python
	```

5. Inspect results under `results/<project_id>/`:

	- `metadata/run_parameters.json` records inputs, QC stats, and (if applicable) the temporary spatial extraction path.
	- `metadata/<project_id>_filter_summary.csv` lists cells before/after each QC threshold.
	- `plots/`, `qc/`, and `tables/` hold spatial/UMAP figures, QC histograms, and clustering tables.
	- `objects/<project_id>_giotto_object.rds` is the serialized Giotto object for downstream use.

Re-run with adjusted thresholds (`--min_genes_per_cell`, `--min_total_expr_per_cell`, `--max_mito_pct`) or a different `--max_cells` value to tune QC and runtime for each dataset.
