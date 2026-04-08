# Quickstart

Use the published GHCR image for the fastest first run. The image is a lean build, so Visium runs should use CSV tissue-position metadata rather than parquet.

1. Pull the published image:

	```bash
	docker pull ghcr.io/nidap-community/giotto-st-pipeline:latest
	```

2. Run a Xenium dataset by mounting inputs at `/data` and results at `/output`:

	```bash
	mkdir -p "$PWD/results/xenium_r1"
	docker run --rm \
	  -v /path/to/xenium:/data:ro \
	  -v "$PWD/results":/output \
	  ghcr.io/nidap-community/giotto-st-pipeline:latest \
	  --stage all \
	  --input_format xenium \
	  --input_dir /data/output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834 \
	  --output_dir /output/xenium_r1 \
	  --project_id XETG00202_R1
	```

3. Run a Visium / Visium HD Spaceranger outs directory. Downsample large runs with `--max_cells`:

	```bash
	docker run --rm \
	  -v /path/to/visium-outs:/data:ro \
	  -v "$PWD/results":/output \
	  ghcr.io/nidap-community/giotto-st-pipeline:latest \
	  --stage all \
	  --input_format visium \
	  --input_dir /data \
	  --output_dir /output/visium_sample \
	  --project_id visium_sample \
	  --max_cells 6000
	```

4. Run an AnnData file:

	```bash
	docker run --rm \
	  -v /path/to/h5ad:/data:ro \
	  -v "$PWD/results":/output \
	  ghcr.io/nidap-community/giotto-st-pipeline:latest \
	  --stage all \
	  --input_format h5ad \
	  --input_path /data/sample.h5ad \
	  --output_dir /output/sample_h5ad \
	  --project_id sample_h5ad \
	  --python_path ~/.local/share/r-miniconda/envs/giotto_env/bin/python
	```

5. Pull the same published image for Apptainer/Singularity when Docker is not available:

	```bash
	singularity pull giotto-st-pipeline.sif docker://ghcr.io/nidap-community/giotto-st-pipeline:latest
	singularity run \
	  --bind /path/to/xenium:/data,/path/to/results:/output \
	  giotto-st-pipeline.sif \
	  --stage all \
	  --input_format xenium \
	  --input_dir /data/output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834 \
	  --output_dir /output/xenium_r1 \
	  --project_id XETG00202_R1
	```

6. Restart from saved workflow checkpoints when needed:

	```bash
	docker run --rm \
	  -v /path/to/results:/output \
	  ghcr.io/nidap-community/giotto-st-pipeline:latest \
	  --stage analyze \
	  --input_object /output/sample_h5ad/objects/sample_h5ad_qc_giotto.rds \
	  --output_dir /output/sample_h5ad \
	  --project_id sample_h5ad
	```

7. Inspect results under `results/<project_id>/`:

	- `metadata/run_parameters.json` records inputs, QC stats, and temporary spatial extraction paths.
	- `metadata/<project_id>_filter_summary.csv` lists cells before and after each QC threshold.
	- `plots/`, `qc/`, and `tables/` hold spatial or UMAP figures, QC histograms, and clustering tables.
	- `objects/<project_id>_ingested_giotto.rds`, `*_qc_giotto.rds`, and `*_analyzed_giotto.rds` are valid restart points.
	- `objects/<project_id>_giotto_object.rds` is the exported final Giotto object for downstream use.

8. Use a local source checkout only when you need to rebuild or debug the workflow itself:

	```bash
	git clone https://github.com/NIDAP-Community/giotto-st-pipeline.git
	cd giotto-st-pipeline
	module load R/4.4.3
	Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
	Rscript -e 'renv::restore(prompt = FALSE)'
	```

Re-run with adjusted thresholds (`--min_genes_per_cell`, `--min_total_expr_per_cell`, `--max_mito_pct`) or a different `--max_cells` value to tune QC and runtime for each dataset.
