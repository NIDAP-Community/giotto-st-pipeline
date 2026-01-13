# File: .github/copilot-instructions.md

You are GitHub Copilot working in this repository to develop a Spatial Transcriptomics pipeline using the R package Giotto, runnable inside a Singularity/Apptainer container.

## 0) First principles (non-negotiable)
- Reproducibility first: any pipeline command shown must be runnable via the container.
- Prefer small, testable increments over large changes.
- Keep changes reversible, especially while debugging.

## 1) Repository bootstrap / location
- This project must live under the GitHub org: https://github.com/NIDAP-Community
- If the remote repo does not exist yet, instruct the user to create it (or request permissions) before suggesting multi-file refactors that assume a remote.
- Default repo name suggestion (unless user overrides): `giotto-st-pipeline` (short, descriptive, searchable).

## 2) Required repository structure (keep standard + predictable)
Maintain (or create) this structure and keep it consistent:

- README.md                  # minimal + runnable
- QUICKSTART.md              # fastest path to first successful run
- CHANGELOG.md               # major milestones only
- container/
  - Singularity.def          # container recipe (primary)
  - build.sh                 # builds .sif deterministically
  - run.sh                   # entrypoint wrapper(s)
  - README.md                # how to build/run container
- scripts/                   # CLI-oriented pipeline scripts (Rscript)
- R/                         # reusable R functions (no side effects)
- configs/                   # config files (YAML/JSON), no secrets
- docs/
  - decision_log.md          # append-only decisions, rationale
  - session_notes.md         # append-only daily notes (>=1x/day)
  - architecture.md          # optional: pipeline overview + data flow
- tests/                     # lightweight checks; prefer testthat if feasible
- examples/                  # tiny example inputs or synthetic generator (small files only)

Do not commit large datasets. Provide documented "data access expectations" instead.

## 3) Decision awareness (MANDATORY)
- Before proposing any major change (architecture, workflow tool, container strategy, data layout), consult `docs/decision_log.md`.
- Do NOT reintroduce previously rejected approaches unless explicitly asked.
- If a suggested change conflicts with a logged decision:
  - flag the conflict explicitly,
  - explain why reconsideration may be warranted,
  - propose an updated decision-log entry rather than silently changing direction.
- When a new decision is made, propose an entry for `docs/decision_log.md`.

## 4) Session notes discipline (MANDATORY)
- `docs/session_notes.md` must be updated at least 1x/day.
- Session notes are append-only: never rewrite history; only append with a new timestamped section.
- At the end of a working session, explicitly ask the user to:
  - `git status`
  - commit with a meaningful message
  - push to GitHub (NIDAP-Community remote)

## 5) Documentation rules (MANDATORY)
### README.md must be current, minimal, and runnable
Keep it short, accurate, and runnable with:
- Setup (container build + prerequisites)
- Data access expectations (what user must provide, where, format)
- How to run the pipeline (single command preferred)
- How to reproduce key outputs (what files are created, where)

### QUICKSTART.md
- A “happy path” that gets a user from zero → first successful run with example/synthetic data.

### CHANGELOG.md
- Add entries only for major milestones (first runnable container, first end-to-end pipeline run, major new dataset support, etc.).

### Removing files
When asked to remove files:
- verify they’re unused,
- search for references across repo,
- if risky, propose a deprecation note (and delay deletion) rather than deleting immediately.

## 6) Container + runtime expectations (MANDATORY)
#### Execution mode (current project decision)
- Current development mode is **local-first** for faster iteration.
- Containerization (Singularity/Apptainer) is a planned milestone once core workflow stabilizes.
- Copilot must consult `docs/decision_log.md` and follow this decision.

### Local reproducibility guardrails
- Maintain an `renv.lock` (or equivalent) and keep it updated with meaningful changes.
- All pipeline entrypoints must remain `Rscript`-driven and non-interactive.
- Avoid OS-specific assumptions; document required system libraries.

### Containerization readiness
- Even during local-first development, design paths, configs, and IO so they can be containerized later
  (no hardcoded user paths; all inputs/outputs configurable).
- Document any system dependencies that must be present in the container later.
- When asked to add container-related instructions, refer to `.github/instructions/container.instructions.md`.

### Container execution expectations
When working on container-related code or instructions:
- All execution examples must use `singularity` or `apptainer` consistently.
- Assume HPC-like constraints unless user says otherwise:
  - limited/no outbound internet during runtime
  - read-only container filesystem
  - bind-mounted working directory
- Prefer containerized R package installation (baked into image) over runtime installs.
- Keep container builds deterministic: pin versions where practical, record them in docs.

- When the container exists, all execution examples must work via Singularity/Apptainer.
- Until then, documentation should present local execution as canonical and container execution as planned.


## 7) Pipeline expectations (Giotto)
When implementing Giotto-based workflow, design scripts/functions around these stages:
- Ingest: define supported input types and required folder layout (document it).
- QC: basic spot/cell filtering, gene filtering; save QC summaries.
- Normalization + feature selection: document defaults and rationale.
- Dimensionality reduction + clustering: parameters configurable.
- Spatial analyses/visualization: at least one reproducible spatial plot output.
- Markers + summary outputs: write tables to disk, not just in-memory objects.
- Optional: cell-type annotation hooks (clearly marked optional).

Always write outputs to a predictable `results/` tree (or configurable output dir) and avoid interactive-only workflows.

## 8) Debugging protocol (MANDATORY order)
When debugging, follow this order:
a) Reproduce the error with the smallest example possible.
b) Add informative errors/warnings where needed (`rlang::abort()`, `cli::cli_abort()`).
c) Add targeted logging via `cli` or `message()` (avoid noisy prints).
d) Fix the root cause; do not patch symptoms.
e) Add a test or cheap check to prevent regression.

If you enter “autonomous/automatic mode”:
- Make incremental changes and keep them reversible.
- Do not perform large refactors while debugging unless asked.
- Summarize what changed and why in a brief changelog-style comment.

## 9) Style + quality for R code
- Prefer functions in `R/` with explicit inputs/outputs; scripts in `scripts/` orchestrate.
- Validate inputs early; fail fast with actionable messages.
- Use consistent paths (prefer `here` or explicit project-relative paths); never hardcode user-specific paths.
- Avoid global state and hidden side effects; write explicit outputs.

## 10) What to do before claiming “done”
Before declaring something complete, ensure:
- The container builds successfully.
- A minimal example run works end-to-end in the container.
- README/Quickstart reflect the exact commands that work.
- Decision log and session notes updated appropriately.

## 11) R script availability
- All R scripts must be runnable via `Rscript` from the command line.
- Avoid interactive-only scripts; if interactivity is needed, provide a non-interactive alternative.
- At a session start time, load R by using `module load R` 
- When terminal is unresponsive to `Rscript`, suggest checking module load status or using full path to Rscript binary.

## 12) Execution model
## Container-as-a-Function (MANDATORY DESIGN INVARIANT)

This project follows a Container-as-a-Function (CAF) model:

- The container exposes a single, stable entrypoint.
- Execution is non-interactive and script-driven.
- Inputs (configs, data paths, flags) are passed explicitly.
- Outputs are written to a user-specified output directory.
- No hidden state or container-internal assumptions are allowed.

Copilot must:
- Preserve entrypoint semantics even during local-first development.
- Avoid designs that require users to “exec into” the container to run analysis.
- Ensure all scripts can be called programmatically via `Rscript`.
- Treat the container as a callable function, not a development shell.

This invariant originates from the `multi-gene-correlations` project and must not be violated without an explicit decision-log update.

## 13) Session start protocol (preferred)
At the start of a coding session (or before major work), Copilot should:
- Review `docs/decision_log.md`, `docs/session_notes.md`, and `docs/architecture.md`
- Summarize relevant constraints + current direction in <10 bullets
- Propose a small, reversible plan
- Wait for confirmation before large changes
- Ensure the user has loaded R via `module load R` or equivalent
- If the terminal is unresponsive to `Rscript`, suggest checking module load status or using full path to Rscript binary.

## 14) When to commit
- Commit at least once per day during active development.
- Commit after completing a logical unit of work (feature, fix, refactor).
- Write meaningful commit messages summarizing changes.
- Make small, focused commits rather than large, monolithic ones.
- Make updates to documents (README, decision log, session notes) part of regular commits.
- Push to GitHub at least once per day to back up work and enable collaboration.
- Do not push autonomous changes without user review.

---



