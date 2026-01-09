---
applyTo: "scripts/run_all.R"
---

Implement `scripts/run_all.R` as the single canonical entrypoint following the **multi-gene-correlations** pattern:

## Entrypoint contract (Container-as-a-Function)
- Must be runnable via: `Rscript scripts/run_all.R [flags]`
- Must be non-interactive; never require user prompts.
- All inputs must be passed explicitly as CLI flags and/or a config file.
- All outputs must be written under `--output_dir` (default `/output`), with a stable folder layout.
- Must create output directories if missing.
- Must write a machine-readable run record (e.g. JSON) to `<output_dir>/metadata/run_parameters.json`.

## Input formats (first-class, explicit)

`run_all.R` must treat the following input formats as **first-class and distinct**:

- `matrix`  
  Explicit expression matrix + spatial coordinates (+ optional metadata).

- `visium`  
  10x Visium directory-based input using Giotto Visium helpers.

- `xenium`  
  10x Xenium output directories containing:
  - `cell_feature_matrix.h5`
  - `cells.csv` or `cells.csv.gz`
  Xenium ingest must NOT be forced into Visium assumptions.

- `auto_dir`  
  A standardized, project-defined “general ST” directory layout (documented separately).

Copilot must:
- Implement `xenium` as its own ingest branch.
- Use `R/Giotto.R` as the reference implementation for Xenium logic.
- Avoid reusing Visium-specific assumptions or helpers for Xenium.

## CLI behavior
- Use a robust CLI parser (`optparse` or `argparse`).
- Accept BOTH:
  1) explicit file paths (expression matrix + spatial coords + metadata)
  2) a standardized input directory layout (general ST support)
- Allow `--config` (YAML/JSON) to specify defaults, with CLI flags overriding config.
- Validate early; fail fast using `cli::cli_abort()` / `rlang::abort()` with actionable messages.
- Use targeted logging via `cli` (no noisy printing).

## General ST format support (Giotto-ready)
- Implement an ingest layer that can harmonize multiple ST formats into a Giotto object:
  - expression matrix (genes x spots/cells) + spatial coordinates + optional metadata
  - Visium should be supported via a dedicated branch, but do NOT assume Visium-only.
- After ingest, downstream steps must be data-source agnostic.

## Output conventions
Create stable subfolders under `--output_dir`:
- `metadata/` (run parameters, session info)
- `objects/` (Giotto object as RDS)
- `qc/` (QC tables/plots)
- `plots/` (spatial plots)
- `tables/` (markers, cluster summaries)
- `logs/` (optional structured log)

## Xenium-specific expectations

When `--input_format xenium` is used:

- `--input_dir` may point to:
  - a single Xenium `output-*` directory, OR
  - a parent directory containing multiple `output-*` directories
- The script must:
  - locate `cell_feature_matrix.h5`
  - locate `cells.csv` or `cells.csv.gz`
  - validate cell/spot ID consistency (as in `R/Giotto.R`)
- Python integration (if required by Giotto) must be parameterized via `--python_path`
  and must not be hardcoded.

## Container-as-a-Function reminder

All input formats, including `xenium`, must obey the Container-as-a-Function contract:
- explicit inputs via flags/config
- explicit outputs under `--output_dir`
- no interactive execution
- no hidden container-internal state

## Local vs container
- The same script must work locally and in-container.
- Never hardcode absolute paths other than default `/output`.
- Assume container runs with bind mounts `/data` and `/output` as common conventions, but do not require them.

## Debugging protocol
When addressing failures, follow the repo debugging order:
(1) minimal repro, (2) informative errors, (3) targeted logging, (4) fix root cause, (5) add a cheap regression check.

## Preflight requirement for run_all changes
Before implementing significant changes to `scripts/run_all.R`, summarize the relevant constraints from:
- `docs/decision_log.md`
- `docs/architecture.md`
Then proceed with minimal reversible changes.


Consult `docs/decision_log.md` before major changes.
