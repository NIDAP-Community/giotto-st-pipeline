# File: .github/instructions/r.instructions.md
---
applyTo: "{R/**/*.R,scripts/**/*.R,**/*.R}"
---

For R in this repo:
- Use `cli` and/or `rlang` for clear errors and targeted logging (no noisy prints).
- Put reusable, testable functions in `R/`; keep `scripts/` as thin orchestrators.
- Validate inputs at the top of functions; fail fast with actionable messages.
- Avoid interactive-only flows; scripts must run via `Rscript` inside the container.
- If you change defaults/parameters, update README/Quickstart and (if major) add a CHANGELOG entry.
