# Decision Log (append-only)

## 2026-01-08 — Initial project conventions
**Decision:** Use `.github/copilot-instructions.md` + optional `.github/instructions/*.instructions.md` for Copilot guidance.  
**Rationale:** Keeps repo-wide rules consistent while allowing file-type-specific guidance.  
**Consequences:** Avoid conflicting instructions across files; resolve conflicts explicitly if they arise.

## 2026-01-08 — Container-as-a-Function execution model
**Decision:** The pipeline container will behave as a pure function: a single, stable entrypoint consumes explicit inputs (arguments + mounted files) and produces explicit outputs, with no hidden state.

**Rationale:** This model maximizes reproducibility, HPC portability, testability, and clarity of execution semantics.

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

## 2026-01-13 — Visium spatial archives + per-filter QC summaries
**Decision:** Treat Spaceranger `*spatial*.zip` bundles as first-class Visium inputs and persist per-threshold QC statistics to a CSV alongside run metadata.

**Rationale:** Many Cell Ranger deliveries package spatial assets as ZIP archives; auto-extracting them preserves the Container-as-a-Function contract. Recording before/after counts for each QC filter improves traceability when tuning thresholds.

**Consequences:**
- `ingest_visium_hd()` now unpacks spatial archives to a temp directory and records provenance in `run_parameters.json`.
- Every run writes `<project_id>_filter_summary.csv` under `metadata/`, enabling downstream reports without parsing JSON.
- Documentation must note archive handling and the additional metadata artefact.

## 2026-04-02 — Align local and container workflow with standard CAF conventions
**Decision:** Use explicit `/data` and `/output` binds, a thin container entrypoint, a helper script for Apptainer/Singularity runs, and explicit `renv` activation for local interactive work instead of repo-level auto-activation.

**Rationale:** This keeps execution semantics consistent across NIDAP CAF repositories, avoids `.Rprofile` side effects during Docker builds, and makes local, Docker, and Apptainer usage easier to document and test consistently.

**Consequences:**
- `.Rprofile` no longer auto-sources `renv/activate.R`; local users must activate `renv` explicitly when working interactively.
- Container docs should describe the `/data` and `/output` bind layout as the canonical interface.
- The repo should provide a host-side Apptainer/Singularity helper to standardize binds and runtime environment.

## 2026-04-02 — Keep one public entrypoint, add stage-based workflow restarts
**Decision:** Preserve a single canonical public entrypoint (`scripts/run_all.R`) while exposing workflow stages through a `--stage` selector and checkpointed Giotto objects.

**Rationale:** The ST pipeline is a workflow with multiple steps, but fragmenting the public container interface into multiple entrypoint scripts would complicate reproducibility, documentation, and HPC execution. Stage-based execution preserves the CAF contract while allowing restarts and checkpointed runs.

**Consequences:**
- `scripts/run_all.R` accepts stage values such as `validate`, `ingest`, `qc`, `analyze`, `export`, and `all`.
- Intermediate Giotto objects are persisted under `objects/` for restartable execution.
- Architecture and user-facing docs must describe stage-aware execution while keeping a single public container entrypoint.

## 2026-04-02 — Publish a lean GHCR image first; document enterprise CA rebuilds for maintainers
**Decision:** Treat the first published GHCR image as a lean build that excludes the optional `arrow` package, and document enterprise CA handling only as a rebuild-time maintainer concern.

**Rationale:** The lean image is more reliable to build and publish, especially in constrained Docker Desktop environments and behind enterprise TLS interception. End users pulling a published image should not need institution-specific CA setup, while maintainers rebuilding from source still need a documented path for internal root CA injection.

**Consequences:**
- The default published image assumes non-parquet Visium metadata at runtime.
- Public quickstart docs should avoid institution-specific certificate instructions for ordinary image consumers.
- Source-build docs should explain BuildKit secret-based CA injection and explicitly forbid committing internal CA files into the repository.

## 2026-04-02 — Add OCI labels before GHCR publication
**Decision:** Add OCI image metadata labels to the Dockerfile before the first GHCR publication.

**Rationale:** GHCR package pages and downstream tooling display OCI metadata directly. Adding the labels before publication improves discoverability and keeps the published image self-describing.

**Consequences:**
- The Dockerfile must carry stable `org.opencontainers.image.*` labels for title, description, source, URL, and documentation.
- Published tags can be inspected without consulting local build context.

## 2026-04-02 — Use one immutable GHCR tag plus floating aliases
**Decision:** Publish each lean image with an immutable `sha-<git-sha>` tag plus floating `lean` and `latest` tags.

**Rationale:** The immutable tag preserves traceability to an exact commit, while the floating tags provide a stable pull target for end users and docs.

**Consequences:**
- GHCR publication should push at least three tags per lean release: `sha-<git-sha>`, `lean`, and `latest`.
- Docs should treat `latest` as the standard pull target while maintainers can pin exact builds with `sha-<git-sha>`.

## 2026-04-08 — Treat the published GHCR image as the default user entrypoint
**Decision:** Make the published GHCR image the default documented entrypoint for routine use, while keeping local source rebuilds as a maintainer workflow.

**Rationale:** The published image now pulls and runs successfully, so end users should not need a source checkout or local `renv` restore just to execute the workflow.

**Consequences:**
- README and QUICKSTART should lead with `docker pull` and `singularity pull` examples.
- Source rebuild instructions should remain in `container/README.md` for maintainers.

## 2026-04-08 — Use GitHub-hosted builds when local rebuilds are blocked
**Decision:** Provide a manual GitHub Actions workflow for lean GHCR publication so maintainers can rebuild and publish the image on GitHub-hosted runners when local Docker builds are blocked by enterprise network behavior.

**Rationale:** Local rebuilds in this environment are currently vulnerable to Ubuntu mirror signature and TLS issues. A GitHub-hosted workflow gives the project a reproducible publication path that is independent of the local workstation network.

**Consequences:**
- The repository should keep a manual GHCR publication workflow under `.github/workflows/`.
- Maintainers can use local Docker builds for debugging, but publication should fall back to GitHub Actions when local rebuilds are unreliable.

## 2026-06-13 — Adopt explicit semantic versioning via VERSION file
**Decision:** Introduce a `VERSION` file at repo root as the single source of truth for pipeline version. The runtime reads it at startup, the Dockerfile bakes it into an OCI label and env var, and git tags must match VERSION for releases.

**Rationale:** Previously, pipeline identity was tracked only via git SHA, which is opaque to end users, silently degrades to "unknown" in non-git contexts (containers, tarballs), and provides no user-facing indicator of compatibility or release status. A semantic version gives users, publications, and CI a stable human-readable identifier.

**Alternatives considered:**
- DESCRIPTION file (R package style): heavier, implies package semantics this project doesn't need.
- Git-tag-only: fails in non-git contexts (container runtime, downloaded archives).
- Env-var-only: requires build-time coordination and has no source-tree presence.

**Consequences:**
- Release workflow must bump `VERSION`, tag the commit, and rebuild the container.
- `provenance.json` now includes `pipeline.version`.
- Runtime warns if version resolves to "unknown" (missing VERSION file + missing env var).
- `--version` flag allows quick identification without a full run.

## 2026-06-13 — Version-tagged container images; remove local publish.sh
**Decision:** Publish container images with a semantic version tag (e.g., `0.1.0`) in addition to `sha-<commit>` and `latest`. Remove the local `container/publish.sh` script since the GitHub Actions workflow is the canonical publish path.

**Rationale:** Users pulling containers could only use opaque SHA tags or the unstable `latest`. A version tag makes the pull command readable and citeable (e.g., `docker pull ...giotto-st-pipeline:0.1.0`). The local `publish.sh` duplicated functionality already covered by the GH Actions workflow and added confusion about which publish path to use.

**Consequences:**
- GHCR images now carry three tags per release: `<version>`, `sha-<commit>`, `latest`.
- README/QUICKSTART container sections reference version tags as the primary pull method.
- Digest pinning is documented as an optional advanced path for bit-level reproducibility.
- `container/publish.sh` is removed; maintainers use the GH Actions workflow or `docker push` directly.

