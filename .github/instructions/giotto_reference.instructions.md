---
applyTo: "{scripts/run_all.R,R/**/*.R}"
---

Use `R/Giotto.R` as the reference implementation for Xenium-style ingest + baseline Giotto workflow.
When implementing `scripts/run_all.R` (and helper modules), preserve the working logic from `R/Giotto.R`
while adapting it to:
- Container-as-a-Function entrypoint semantics (explicit inputs, explicit outputs, no hidden state)
- Local-first development (but container-ready defaults like --output_dir=/output)
- General ST support via multiple ingest branches (Xenium as one branch, not the only one)

Do NOT copy hardcoded local paths (e.g. `/data/STAG/data`, `~/.local/.../python`).
Instead, parameterize them via CLI flags and/or config.
