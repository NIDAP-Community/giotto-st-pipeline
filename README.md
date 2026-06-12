# Giotto Spatial Transcriptomics Pipeline

Container-ready R workflow for first-pass spatial transcriptomics analysis with Giotto. It ingests Xenium, Visium/Visium HD, and AnnData inputs, then exports reproducible QC plots, embeddings, cluster tables, Giotto objects, and run metadata.

The workflow is designed around a simple container contract:

- mount input data at `/data`
- mount a writable results directory at `/output`
- run one stable entrypoint with explicit command-line options

For source development, rebuilds, dependency maintenance, and Giotto version notes, see [docs/developer_readme.md](docs/developer_readme.md).

## When To Use This Workflow

Use this workflow for standardized screening across spatial transcriptomics datasets. It is meant to answer practical first-pass questions: does the dataset ingest cleanly, what do the QC distributions look like, are there clear UMAP or spatial cluster patterns, and which samples need deeper follow-up?

This workflow is not intended to replace exploratory Giotto analysis. For custom plotting, project-specific interpretation, method development, image-aware review, or newer Giotto Suite features, working directly in Giotto or Giotto Suite is usually the better next step.

## Example Outputs

Visium mouse brain example (`--max_cells 6000`):

![Spatial clusters](docs/images/visium_mouse_brain_readme_spatial.png)
![UMAP embedding](docs/images/visium_mouse_brain_readme_umap.png)
![Genes per cell histogram](docs/images/visium_mouse_brain_readme_nr_genes_hist.png)
![Total expression histogram](docs/images/visium_mouse_brain_readme_total_expr_hist.png)
![Genes vs expression scatter](docs/images/visium_mouse_brain_readme_genes_vs_expr.png)

These figures are exported automatically beneath `results/<project_id>/plots/` and `results/<project_id>/qc/` for every pipeline run.

## Supported Inputs

- Xenium output directories with `cell_feature_matrix.h5` and cells metadata
- Visium and Visium HD Spaceranger `outs/` directories
- AnnData `.h5ad` files

The default published container image is intended to restore the full locked runtime, including optional packages such as `arrow`. For smaller local debugging builds, maintainers may publish or build an explicit lean variant with reduced Visium parquet support.

## What The Pipeline Produces

- QC metrics, filter summaries, and QC plots
- UMAP and spatial cluster plots
- Leiden cluster assignments
- Restartable Giotto object checkpoints after ingest, QC, and analysis stages
- A final exported Giotto object
- `metadata/run_parameters.json` and `metadata/session_info.txt` for provenance

## Quick Start

The recommended path is the published GHCR image.

```bash
docker pull ghcr.io/nidap-community/giotto-st-pipeline:latest
```

For reproducible research, pin the image used in a run. The `latest` tag is
convenient for trying the workflow, but it is intentionally movable. Published
releases include an immutable source tag:

```bash
docker pull ghcr.io/nidap-community/giotto-st-pipeline:sha-<full-git-sha>
```

For the strongest pin, use the registry digest reported by `docker pull` or
`docker inspect`:

```bash
docker pull ghcr.io/nidap-community/giotto-st-pipeline@sha256:<image-digest>
```

Run a Xenium dataset:

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

Run a Visium or Visium HD Spaceranger `outs/` directory:

```bash
mkdir -p "$PWD/results/visium_sample"
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

Run an AnnData `.h5ad` file:

```bash
mkdir -p "$PWD/results/sample_h5ad"
docker run --rm \
	-v /path/to/h5ad:/data:ro \
	-v "$PWD/results":/output \
	ghcr.io/nidap-community/giotto-st-pipeline:latest \
	--stage all \
	--input_format h5ad \
	--input_path /data/sample123.h5ad \
	--output_dir /output/sample_h5ad \
	--project_id sample_h5ad \
	--python_path ~/.local/share/r-miniconda/envs/giotto_env/bin/python
```

For more examples, including restartable stage runs, see [QUICKSTART.md](QUICKSTART.md).

## HPC Apptainer/Singularity

On HPC systems where Docker is not available, pull the published image into a `.sif` file once:

```bash
module load singularity
singularity pull giotto-st-pipeline.sif docker://ghcr.io/nidap-community/giotto-st-pipeline:latest
```

For a fixed image, replace `:latest` with either `:sha-<full-git-sha>` or
`@sha256:<image-digest>`.

Some clusters provide `apptainer` instead of `singularity`; the commands are otherwise equivalent:

```bash
module load apptainer
apptainer pull giotto-st-pipeline.sif docker://ghcr.io/nidap-community/giotto-st-pipeline:latest
```

Run analysis jobs on a compute node, binding input data read-only at `/data` and a writable results directory at `/output`:

```bash
mkdir -p /path/to/results/visium_sample
singularity run --cleanenv \
	--bind /path/to/spaceranger/outs:/data:ro \
	--bind /path/to/results:/output \
	giotto-st-pipeline.sif \
	--stage all \
	--input_format visium \
	--input_dir /data \
	--output_dir /output/visium_sample \
	--project_id visium_sample \
	--max_cells 6000 \
	--spatial_point_size 2.25
```

For Xenium, bind the parent Xenium output directory and set `--input_format xenium`. For `.h5ad`, bind the folder containing the file and use `--input_format h5ad --input_path /data/sample.h5ad`.

## Common Options

| Flag | Description |
| --- | --- |
| `--stage` | Workflow stage: `all`, `validate`, `ingest`, `qc`, `analyze`, or `export`. Default is `all`. |
| `--input_format` | Choose `xenium`, `visium`, or `h5ad`; default `auto` infers from directory structure or file extension. |
| `--input_dir` | Input directory for Xenium or Visium/Visium HD runs. |
| `--input_path` | Direct path to a single-file input, currently `.h5ad`. |
| `--input_object` | Existing Giotto object RDS used as input for restartable `qc`, `analyze`, or `export` stages. |
| `--output_dir` | Directory where results are written. |
| `--project_id` | Short identifier used to prefix plot and table artefacts. Defaults to folder name. |
| `--max_cells` | Randomly subsample cells/spots before analysis; useful for large Visium HD runs. |
| `--min_genes_per_cell` | Drop cells whose detected genes fall below this threshold. |
| `--min_total_expr_per_cell` | Drop cells whose total expression counts fall below this threshold. |
| `--max_mito_pct` | Drop cells whose mitochondrial expression fraction exceeds this percentage. |
| `--mito_gene_prefixes` | Comma-separated mitochondrial gene prefixes; default is `MT-`; pass `none` to disable. |
| `--pca_dims` | Maximum number of PCA dimensions used for UMAP and clustering. Default is `10`. |
| `--neighbor_k` | Nearest-neighbor graph k before clustering. Default is `20`; capped at cells - 1. |
| `--cluster_resolution` | Leiden clustering resolution. Higher values usually produce more clusters. |
| `--spatial_point_size` | Dot size for the exported spatial map. |
| `--umap_point_size` | Dot size for the exported UMAP plot. |
| `--python_path` | Optional Python binary for Giotto and AnnData support. |
| `--cores` | Number of CPU cores to dedicate to Giotto. Default is `4`. |
| `--seed` | Random seed applied before dimensionality reduction and clustering. Default is `1`. |
| `--dry_run` | Validate inputs and exit before running the workflow. |

## Iterating On Results

First-pass clustering may be too coarse or too granular. Re-run with adjusted `--cluster_resolution`, `--pca_dims`, or `--neighbor_k` to tune the preliminary view before deciding whether a dataset needs deeper interactive analysis.

To rerun only plotting and table export from an existing analyzed object, use `--stage export` with `--input_object`. This is useful when adjusting plot settings such as `--spatial_point_size`, legend size, or axis text without repeating ingest, QC, and clustering.

```bash
docker run --rm \
	-v /path/to/results:/output \
	ghcr.io/nidap-community/giotto-st-pipeline:latest \
	--stage export \
	--input_object /output/visium_sample/objects/visium_sample_analyzed_giotto.rds \
	--output_dir /output/visium_sample \
	--project_id visium_sample \
	--spatial_point_size 2.25 \
	--spatial_legend_text 12 \
	--spatial_axis_title 12
```

## Output Layout

```text
results/<project_id>/
├── metadata/
│   ├── run_parameters.json
│   ├── session_info.txt
│   └── <project_id>_filter_summary.csv
├── objects/
│   ├── <project_id>_ingested_giotto.rds
│   ├── <project_id>_qc_giotto.rds
│   ├── <project_id>_analyzed_giotto.rds
│   └── <project_id>_giotto_object.rds
├── qc/
│   ├── <project_id>_qc_metrics.csv
│   ├── <project_id>_qc_summary.txt
│   ├── <project_id>_nr_genes_hist.png
│   ├── <project_id>_total_expr_hist.png
│   └── <project_id>_genes_vs_expr.png
├── plots/
│   ├── <project_id>_spatial.png
│   └── <project_id>_umap.png
└── tables/
    └── clusters.csv
```

When executed inside the published container, the same directory layout is written under the mounted `--output_dir`.

## Notes And Limitations

- Visium runs generated by Spaceranger sometimes package spatial assets inside `*spatial*.zip`; the pipeline extracts these archives automatically during ingest.
- AnnData ingest relies on Python packages `anndata` and `scipy` being available to the Giotto reticulate environment.
- For very high-resolution Visium HD runs, use `--max_cells` or schedule the job on a compute node to avoid memory failures on login nodes.
- Matrix-format ingest remains TODO.

## Additional Documentation

- [QUICKSTART.md](QUICKSTART.md): concise runnable examples
- [docs/data_description.md](docs/data_description.md): input and output data expectations
- [docs/architecture.md](docs/architecture.md): architecture and execution model
- [container/README.md](container/README.md): container build, run, and publication details
- [docs/developer_readme.md](docs/developer_readme.md): local source workflow, rebuilds, dependency maintenance, and roadmap
