# Configuration Templates

YAML or JSON config files are the recommended way to run this pipeline. The
same config file can be used from:

- local command line R with `Rscript scripts/run_all.R --config ...`
- R Markdown notebooks
- Singularity/Apptainer containers
- Docker containers

CLI flags always take precedence over config values, so a saved config can be
reused while changing one option at the command line.

## Templates

```text
configs/
├── xenium_template.yaml
├── visium_template.yaml
├── h5ad_template.yaml
└── README.md
```

Template files:

- [xenium_template.yaml](xenium_template.yaml)
- [visium_template.yaml](visium_template.yaml)
- [h5ad_template.yaml](h5ad_template.yaml)

Copy the closest template, rename it for your dataset, and replace the example
paths before running the workflow. The templates list every common flag so users
can see what is configurable. Use `none` for optional values you want to leave
off. Keep sensitive paths or PHI out of version control; commit generic
templates rather than study-specific paths.

## Parameter Reference

Use `input_dir` for Visium and Xenium runs. Use `input_path` for h5ad runs.

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_format` | string | One of: `visium`, `xenium`, `h5ad`, `matrix`, `auto`. Auto-detects from directory contents if set to `auto`. |
| `input_dir` | path | Directory containing Visium or Xenium data, such as a Spaceranger `outs/` directory or Xenium output folder. |
| `input_path` | path | Direct path to an `.h5ad` file. Use this instead of `input_dir` for h5ad runs. |
| `output_dir` | path | Where results are written. Use `/output` for container runs. |
| `project_id` | string | Short name for output files (e.g., `sample_A`). No spaces or special characters. |

### Optional Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cores` | integer | 4 | Number of CPU cores for parallel steps. |
| `seed` | integer | 1 | Random seed for reproducibility. |
| `max_cells` | integer or `none` | `none` | Randomly downsample to this many cells before analysis. Useful for large Visium HD. |
| `pca_dims` | integer | 10 | Number of PCA dimensions for UMAP and clustering. |
| `neighbor_k` | integer | 20 | k for nearest-neighbor graph (auto-capped at cells - 1). |
| `cluster_resolution` | float | 0.4 | Leiden clustering resolution. Higher = more clusters. |

### Optional QC Thresholds

Use `none` or omit the field to skip a threshold entirely. The templates show
the fields explicitly as `none` so users know the option exists without turning
on filtering by accident.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_genes_per_cell` | integer or `none` | `none` | Drop cells detecting fewer genes than this. |
| `min_total_expr_per_cell` | integer or `none` | `none` | Drop cells with total UMI counts below this. |
| `max_mito_pct` | float or `none` | `none` | Drop cells with mitochondrial fraction above this %. |
| `mito_gene_prefixes` | string | `MT-` | Comma-separated prefixes for mito genes. Use `none` to disable. |

### Plot Appearance

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `spatial_point_size` | float | 2.25 | Dot size on spatial map. |
| `umap_point_size` | float | 1.5 | Dot size on UMAP. |
| `spatial_legend_text` | float | 12 | Legend text size (pt). |
| `spatial_legend_symbol_size` | float | 1.4 | Legend symbol scale factor. |
| `spatial_axis_text` | float | 12 | Axis tick label size (pt). |
| `spatial_axis_title` | float | 12 | Axis title size (pt). |

### Advanced / Workflow Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `stage` | string | `all` | Run only a specific stage: `validate`, `ingest`, `qc`, `analyze`, `export`, or `all`. |
| `input_object` | path or `none` | `none` | Path to a saved Giotto `.rds` for restarting from a checkpoint. |
| `python_path` | path or `none` | `none` | Explicit Python binary for Giotto's reticulate integration. |
| `dry_run` | boolean | false | Validate config and inputs, then exit without running analysis. |
| `verbose` | boolean | false | Enable verbose logging output. |

## Validating Your Config

Before a long run, validate that paths exist and parameters parse correctly:

```bash
Rscript scripts/validate_config.R configs/my_run.yaml
```

This checks input paths, parameter types, and data format detection without
running the full analysis.

## Local CLI

```bash
Rscript scripts/run_all.R --config configs/my_visium_run.yaml
```

Override a single value without editing the config:

```bash
Rscript scripts/run_all.R \
	--config configs/my_visium_run.yaml \
	--cluster_resolution 0.8
```

## Container Runs

For containers, paths in the config must match the mounted container paths. The
standard convention is:

- `/data`: read-only input data and optional config files
- `/output`: writable analysis output

For example, a Visium config used inside a container might set:

```yaml
input_format: visium
input_dir: /data
output_dir: /output/visium_sample
project_id: visium_sample
```

## Publication Provenance

Every run writes:

- `metadata/provenance.json`
- `metadata/run_parameters.json`
- `metadata/session_info.txt`

Archive those files with the results used for publication.
