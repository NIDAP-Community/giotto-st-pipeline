## Data Folder Directives

The test data root is `/Users/maggiec/GitHub/Maggie/NIDAP/SpatialData`.

Treat this folder as local fixture data for development and testing, not as authoritative production data.

Folder layout:
- `xenium/`
  - Contains a lightweight Xenium fixture folder:
    - `output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834/`
  - This fixture includes representative Xenium metadata and matrix files:
    - `experiment.xenium`
    - `metrics_summary.csv`
    - `gene_panel.json`
    - `cells.parquet`
    - `cell_feature_matrix.h5`
    - `cell_feature_matrix/`
    - `analysis_summary.html`
  - This is a trimmed subset, not the full original Xenium output.

- `visium/`
  - Contains one sample-style fixture at:
    - `sample123/outs/spatial/`
  - This folder holds Visium spatial assets only, such as images, scale factors, and tissue-position files.
  - It is structured to resemble a standard Visium sample layout for testing.

- `h5ad/`
  - Contains:
    - `sample123.h5ad` — larger source fixture
    - `sample123_small.h5ad` — small test fixture
  - Prefer `sample123_small.h5ad` for routine tests and examples.
  - `sample123_small.h5ad` has shape `(500, 1000)`.
  - `sample123.h5ad` has shape `(22260, 17943)`.

- `results/`
  - Output directory for generated test artifacts, intermediate files, or pipeline results.
  - Safe place for derived outputs during local development.

Behavioral guidance:
- Prefer the smallest fixture that satisfies the task.
- Do not assume full raw vendor outputs are present.
- Do not overwrite source fixtures in `xenium/`, `visium/`, or `h5ad/` unless explicitly requested.
- Write derived files to `results/` whenever possible.
- Keep tests lightweight and deterministic by using the trimmed fixtures.
