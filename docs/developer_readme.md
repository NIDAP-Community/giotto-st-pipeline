# Developer README

This document collects source-checkout, rebuild, dependency, and maintenance guidance for the Giotto Spatial Transcriptomics Pipeline. End users who only need to run the published image should start with the main [README](../README.md) and [QUICKSTART](../QUICKSTART.md).

## Repository Map

- `scripts/run_all.R`: canonical CLI entrypoint and workflow controller
- `R/`: reusable ingest, QC, Giotto analysis, and export modules
- `container/`: Docker, Apptainer/Singularity, and GHCR publication helpers
- `renv.lock`: pinned R package environment
- `docs/`: architecture notes, decision logs, data descriptions, and development notes
- `configs/`: example and future configuration files

The pipeline follows a Container-as-a-Function model: inputs are mounted under `/data`, outputs are written under `/output`, and all runtime choices are explicit command-line arguments.

## Local R Workflow

`renv` is available, but it is not auto-activated by `.Rprofile`.

```bash
module load R/4.4.3
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore(prompt = FALSE)'
```

For interactive local work, activate the project library explicitly:

```r
source("renv/activate.R")
```

For non-interactive runs, call `Rscript scripts/run_all.R ...` from the repo root after restoring `renv`. Because `.Rprofile` does not auto-activate the project library, either run from an R environment where the restored packages are already visible, or pass the restored `renv` library through `R_LIBS_USER`:

```bash
PROJECT_R_LIB="$(Rscript -e 'source("renv/activate.R"); cat(.libPaths()[1])')"
R_LIBS_USER="$PROJECT_R_LIB" Rscript scripts/run_all.R --help
```

If Giotto should use a specific Python environment, pass it both to reticulate and to the pipeline:

```bash
PYTHON_BIN="/path/to/giotto/python"
```

The examples below use `RETICULATE_PYTHON="$PYTHON_BIN"` and `--python_path "$PYTHON_BIN"` together so both reticulate and the pipeline point at the same Python binary.

## Source Checkout Runs

Run these examples from the repo root after restoring the R environment. Replace the input and output paths with local paths on the machine where the command is running.

```bash
PROJECT_R_LIB="$(Rscript -e 'source("renv/activate.R"); cat(.libPaths()[1])')"
PYTHON_BIN="/path/to/giotto/python"
```

Xenium:

```bash
mkdir -p results/xenium_r1
RETICULATE_PYTHON="$PYTHON_BIN" \
	R_LIBS_USER="$PROJECT_R_LIB" \
	Rscript scripts/run_all.R \
	--stage all \
	--input_format xenium \
	--input_dir /path/to/xenium/output-XETG00202__sample \
	--output_dir results/xenium_r1 \
	--project_id xenium_r1 \
	--python_path "$PYTHON_BIN"
```

Visium or Visium HD Spaceranger `outs/` directory:

```bash
mkdir -p results/visium_sample
RETICULATE_PYTHON="$PYTHON_BIN" \
	R_LIBS_USER="$PROJECT_R_LIB" \
	Rscript scripts/run_all.R \
	--stage all \
	--input_format visium \
	--input_dir /path/to/spaceranger/outs \
	--output_dir results/visium_sample \
	--project_id visium_sample \
	--max_cells 6000 \
	--spatial_point_size 2.25 \
	--python_path "$PYTHON_BIN"
```

AnnData `.h5ad` file:

```bash
mkdir -p results/sample_h5ad
RETICULATE_PYTHON="$PYTHON_BIN" \
	R_LIBS_USER="$PROJECT_R_LIB" \
	Rscript scripts/run_all.R \
	--stage all \
	--input_format h5ad \
	--input_path /path/to/sample.h5ad \
	--output_dir results/sample_h5ad \
	--project_id sample_h5ad \
	--python_path "$PYTHON_BIN"
```

To rerun only plotting and table export from an existing analyzed object, use `--stage export` with `--input_object`:

```bash
RETICULATE_PYTHON="$PYTHON_BIN" \
	R_LIBS_USER="$PROJECT_R_LIB" \
	Rscript scripts/run_all.R \
	--stage export \
	--input_object results/visium_sample/objects/visium_sample_analyzed_giotto.rds \
	--output_dir results/visium_sample \
	--project_id visium_sample \
	--spatial_point_size 2.25 \
	--spatial_legend_text 12 \
	--spatial_axis_title 12 \
	--python_path "$PYTHON_BIN"
```

For clustering iteration, rerun `--stage all` or `--stage analyze` with updated values for `--pca_dims`, `--neighbor_k`, or `--cluster_resolution`. Higher `--cluster_resolution` values usually produce more clusters; lower values usually produce fewer clusters.

## Container Rebuilds

Restore the R environment locally, then build the image on a workstation with Docker:

```bash
./container/build.sh giotto-st-pipeline:dev
```

On Apple Silicon, pass the target platform explicitly when needed:

```bash
./container/build.sh giotto-st-pipeline:dev --platform linux/amd64
```

Optionally export to a tarball and convert to `.sif` for Apptainer/Singularity:

```bash
docker save giotto-st-pipeline:dev -o giotto-st-pipeline.tar
singularity build giotto-st-pipeline.sif docker-archive://giotto-st-pipeline.tar
```

HPC environments without Docker should rely on the published GHCR image or pre-built `.sif` artefacts generated off-cluster.

For standardized Apptainer/Singularity binds, use `container/run_apptainer.sh` with `DATA_DIR` and `OUTPUT_DIR` environment variables.

If you are only consuming the published GHCR image, you should not need to manage enterprise CA certificates locally. Enterprise CA handling is only relevant when rebuilding the image from source behind a TLS-inspecting proxy.

For GHCR publication details, tag conventions, Docker login requirements, optional fuller builds, and enterprise CA handling, see [container/README.md](../container/README.md).

## Reproducible Environment

- `renv.lock` pins CRAN, Bioconductor, and GitHub package revisions, including Giotto and spatstat suites.
- `container/Dockerfile` pins the Rocker base image by digest and pins the `renv` bootstrap version before restore.
- `container/python-bootstrap-requirements.txt` and `container/python-requirements.txt` pin the Python environment used by reticulate-backed ingest and clustering helpers.
- `.Rprofile` does not auto-activate renv; activate it explicitly with `source("renv/activate.R")` for local interactive sessions.
- System dependencies required for compiled R packages are documented in `container/Dockerfile`.

Use `renv::status()` before committing dependency changes to ensure the lockfile stays current. For strict archival container reproducibility, build against a pinned Ubuntu apt snapshot or frozen internal mirror; a digest-pinned base image does not freeze packages installed later by `apt-get update`.

Published images should always include an immutable `sha-<full-git-sha>` tag.
That source commit contains the exact `renv.lock` restored by the container
build. The container build also stamps OCI labels for the source revision and
`renv.lock` SHA-256; keep those labels aligned with the published tag.

## Giotto Version Stance

As of 2026-06-10, this workflow pins legacy `Giotto 1.1.2` from the original `RubD/Giotto` line in `renv.lock`. The current checked Giotto Suite release is `v4.2.3`, released on 2026-06-02.

Intermediate Suite-era releases exist, including `v3.3.2`, but they are not expected to be low-upheaval drop-in upgrades for this repo. Upgrading from `Giotto 1.1.2` should be treated as a compatibility project because the Suite-era releases changed object structures, accessors, import paths, and package organization.

For standardized preliminary analyses, the pinned `Giotto 1.1.2` environment is acceptable unless a specific newer Giotto Suite capability is needed. Keep the stable workflow pinned for routine screening, and evaluate Giotto Suite separately when deeper or more modern spatial analysis features are required.

Version references:

- Legacy Giotto `1.1.2`: <https://raw.githubusercontent.com/RubD/Giotto/master/DESCRIPTION>
- Current Giotto Suite release history: <https://github.com/giotto-suite/Giotto/releases>
- Current Giotto Suite `DESCRIPTION`: <https://raw.githubusercontent.com/giotto-suite/Giotto/suite/DESCRIPTION>

## Maintenance Notes

- `scripts/run_all.R` is the canonical CLI hub. Keep new user-facing knobs wired through this entrypoint.
- `R/pipeline_basic.R` owns the shared analysis and export behavior.
- `container/Dockerfile` restores from `renv.lock`; avoid maintaining a separate hand-curated R package list.
- Publication uses the GitHub Actions workflow (`.github/workflows/publish-ghcr.yml`); it reads `VERSION` and publishes version, SHA, and `latest` tags.
- The default container build restores the full lockfile. Use `--build-arg RENV_RESTORE_EXCLUDE=arrow` only for an explicit lean/debug variant.
- `results/` is ignored by git; generated analysis outputs should be archived externally when needed.

## Roadmap

- Add matrix-format ingest helper.
- Decide whether to publish a second fuller image variant with parquet-enabled Visium support.
- Automate lightweight tests under `tests/`.
- Document example configs under `configs/`.
