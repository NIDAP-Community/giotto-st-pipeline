---
applyTo: "docs/architecture.md"
---

When editing or extending `docs/architecture.md`:

## Purpose
- Treat this file as a **high-level system map**, not a user manual.
- Explain *how components fit together*, *how data flows*, and *where responsibilities live*.
- Assume the reader is a technical scientist (R/ST/omics-aware) but new to this repo.

## Scope rules
Include:
- Major components (container, scripts, R modules, configs, docs)
- Execution model (entrypoints, orchestration, independence of stages)
- Data flow (inputs → transformations → outputs)
- Supported data *classes* (general ST formats), not specific datasets
- Container/runtime assumptions (HPC, no internet, bind mounts)

Exclude:
- Step-by-step commands (belongs in README/Quickstart)
- Parameter exhaustiveness (belongs in configs docs)
- Daily progress or experiments (belongs in session notes)
- Justification debates (belongs in decision_log.md)

## Decision awareness
- If architectural structure changes, reference the relevant decision in `docs/decision_log.md`.
- Do not contradict logged decisions.
- If an architectural change implies a new decision, propose a decision-log entry explicitly.

## Style
- Prefer diagrams (ASCII is fine) + concise prose.
- Use stable concepts; avoid volatile implementation details.
- Keep the file readable in <10 minutes.

## Maintenance
- Update this file only for **structural changes**, not bug fixes.
- When updated, ensure README/Quickstart remain consistent at a high level.
