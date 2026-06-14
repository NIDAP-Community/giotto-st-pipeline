# Quickstart

The same pipeline can be run from local R, R Markdown, Singularity/Apptainer, or
Docker. For most R users, start with a local checkout and a YAML config file.

## 1. Local Command-Line R

```bash
git clone https://github.com/NIDAP-Community/giotto-st-pipeline.git
cd giotto-st-pipeline

# For publication, pin to a specific release before restoring packages:
# git checkout v0.1.0

Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore(prompt = FALSE)'

export R_LIBS_USER="$(Rscript -e 'source("renv/activate.R"); cat(.libPaths()[1])')"
cp configs/visium_template.yaml configs/my_visium_run.yaml
```

Edit `configs/my_visium_run.yaml` so `input_dir`, `output_dir`, and
`project_id` match your dataset. For h5ad runs, use `input_path` instead of
`input_dir`. Validate before running:

```bash
Rscript scripts/validate_config.R configs/my_visium_run.yaml
```

If validation passes, run the pipeline:

```bash
Rscript scripts/run_all.R --config configs/my_visium_run.yaml
```

Use the same pattern for:

- [configs/xenium_template.yaml](configs/xenium_template.yaml)
- [configs/visium_template.yaml](configs/visium_template.yaml)
- [configs/h5ad_template.yaml](configs/h5ad_template.yaml)

List all command-line flags:

```bash
Rscript scripts/run_all.R --help
```

Override a saved config value without editing the file:

```bash
Rscript scripts/run_all.R \
	--config configs/my_visium_run.yaml \
	--cluster_resolution 0.8
```

## 2. R Markdown

Use the example notebook when you want the command, outputs, and provenance in a
single rendered report.

```bash
Rscript -e 'rmarkdown::render("examples/run_pipeline.Rmd", params = list(config = "configs/my_visium_run.yaml"))'
```

The notebook calls `scripts/run_all.R --config ...`, then displays paths to:

- `metadata/provenance.json`
- `metadata/run_parameters.json`
- `metadata/session_info.txt`
- QC metrics
- cluster table
- spatial and UMAP plots

## 3. Singularity Or Apptainer

On HPC, pull the container image:

```bash
module load singularity

IMAGE_REPO="ghcr.io/nidap-community/giotto-st-pipeline"

# Use "latest" for testing, or a release version for publication:
IMAGE_TAG="0.1.0"

singularity pull giotto-st-pipeline.sif "docker://${IMAGE_REPO}:${IMAGE_TAG}"
```

Create a container-facing config where paths use `/data` and `/output`, for
example:

```bash
cp configs/visium_template.yaml configs/my_visium_container.yaml
```

Edit `configs/my_visium_container.yaml`:

```yaml
stage: all
input_format: visium
input_dir: /data
output_dir: /output/visium_sample
project_id: visium_sample
max_cells: 6000
```

Then run:

```bash
DATA_DIR=/path/to/spaceranger/outs \
OUTPUT_DIR=/path/to/results \
EXTRA_BINDS="$PWD/configs:/configs:ro" \
./container/run_apptainer.sh giotto-st-pipeline.sif \
	--config /configs/my_visium_container.yaml
```

The helper script works with either `apptainer` or `singularity`, mounts
`DATA_DIR` read-only at `/data`, mounts `OUTPUT_DIR` at `/output`, and passes
the remaining flags to `scripts/run_all.R` inside the container.

List all command-line flags from the container:

```bash
./container/run_apptainer.sh giotto-st-pipeline.sif --help
```

## 4. Docker

Use Docker if it is already part of your workstation workflow.

```bash
IMAGE_REPO="ghcr.io/nidap-community/giotto-st-pipeline"

# Use "latest" for testing, or a release version for publication:
IMAGE_TAG="0.1.0"
docker pull "${IMAGE_REPO}:${IMAGE_TAG}"

cp configs/visium_template.yaml configs/my_visium_container.yaml
# Edit configs/my_visium_container.yaml so input_dir is /data and output_dir is /output/...

docker run --rm \
	-v /path/to/spaceranger/outs:/data:ro \
	-v "$PWD/results":/output \
	-v "$PWD/configs":/configs:ro \
	"${IMAGE_REPO}:${IMAGE_TAG}" \
	--config /configs/my_visium_container.yaml
```

Check what's inside the image:

```bash
docker inspect "${IMAGE_REPO}:${IMAGE_TAG}" \
	--format 'version={{ index .Config.Labels "org.opencontainers.image.version" }} commit={{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

List all command-line flags from Docker:

```bash
docker run --rm "${IMAGE_REPO}:${IMAGE_TAG}" --help
```

## Record For Publication

Every run writes publication metadata under `results/<project_id>/metadata/`:

- `provenance.json`
- `run_parameters.json`
- `session_info.txt`

## Inspect Results

Outputs are written under `results/<project_id>/` or the `output_dir` from your
config:

- `metadata/provenance.json` records pipeline, lockfile, runtime, and container provenance.
- `metadata/run_parameters.json` records inputs, options, QC stats, and output paths.
- `metadata/<project_id>_filter_summary.csv` lists cells before and after each QC threshold.
- `plots/`, `qc/`, and `tables/` hold spatial or UMAP figures, QC histograms, and clustering tables.
- `objects/<project_id>_ingested_giotto.rds`, `*_qc_giotto.rds`, and `*_analyzed_giotto.rds` are restart points.
- `objects/<project_id>_giotto_object.rds` is the exported final Giotto object for downstream use.
