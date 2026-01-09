# Decision Log (append-only)

## 2026-01-08 — Initial project conventions
**Decision:** Use `.github/copilot-instructions.md` + optional `.github/instructions/*.instructions.md` for Copilot guidance.  
**Rationale:** Keeps repo-wide rules consistent while allowing file-type-specific guidance.  
**Consequences:** Avoid conflicting instructions across files; resolve conflicts explicitly if they arise.

## 2026-01-08 — Container-as-a-Function execution model
**Decision:** The pipeline container will behave as a pure function: a single, stable entrypoint consumes explicit inputs (arguments + mounted files) and produces explicit outputs, with no hidden state.

**Rationale:** This model proved robust in the `multi-gene-correlations` project. It maximizes reproducibility, HPC portability, testability, and clarity of execution semantics.

**Constraints / Invariants:**
- One canonical entrypoint script (e.g., `run_all.R`)
- No interactive execution required
- All inputs passed explicitly (config files, paths, flags)
- All outputs written to a user-specified output directory
- No reliance on container-internal writable state

**Consequences:**
- Development (local or containerized) must preserve entrypoint semantics
- Containers are callable like functions, not treated as interactive environments
- Pipeline scripts must be composable and side-effect aware

## 2026-01-09 — Adopt renv for dependency capture
**Decision:** Snapshot R dependencies with renv (`renv.lock`) and treat renv restore as the canonical way to reproduce the local environment.

**Rationale:** Locking package versions avoids "works on my machine" drift, supports deterministic Docker/Singularity builds, and mirrors the reproducibility goals from previous CAF projects.

**Consequences:**
- Developers must run `renv::restore()` (or use the container) before executing scripts.
- New R packages require an updated `renv.lock`; large install steps should occur outside the HPC head node when possible.
- Docker/Apptainer images will install dependencies via `renv::restore()` during build time.

