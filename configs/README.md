# Configuration Templates

Use this directory to store runnable configuration bundles for `scripts/run_all.R`. Each file should be
self-documenting (YAML or JSON) and correspond to a specific dataset, QC profile, or container invocation.

## Recommended Structure

```
configs/
├── xenium_default.yaml       # Minimal flags for standard Xenium runs
├── xenium_highstringency.yaml # Example filtering thresholds (TODO)
└── README.md
```

### Example: `xenium_default.yaml`

```yaml
input_format: xenium
input_dir: /data/xenium/output-EXAMPLE_RUN
output_dir: /output/xenium_example
project_id: EXAMPLE_RUN
cores: 4
seed: 1
```

### Usage

```bash
Rscript scripts/run_all.R --config configs/xenium_default.yaml
```

- CLI flags always take precedence over config values.
- Keep sensitive paths (e.g., PHI) out of version control; provide templates or placeholders instead.
- Update this README whenever a new canonical config is added.
