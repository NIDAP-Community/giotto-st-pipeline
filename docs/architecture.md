# Architecture Overview

## Purpose

This document describes the high-level architecture of the **Giotto Spatial Transcriptomics (ST) Pipeline**.
It explains how major components fit together, how data flows through the system, and the execution model
used to ensure reproducibility in HPC environments.

This document is intentionally **implementation-light**. Detailed commands, parameters, and examples live
in the README and Quickstart.

---

## Design Goals

- Support **general Spatial Transcriptomics formats** compatible with Giotto
- Run reproducibly inside a **Singularity/Apptainer container**
- Be modular, script-driven, and suitable for HPC environments
- Avoid interactive-only or notebook-centric workflows
- Make data flow and responsibilities explicit

---

## Supported Data Model (General ST)

The pipeline is designed around **Giotto’s abstract data model**, not a single vendor platform.

Supported input classes include:
- Visium-style spot-based ST
- Other grid- or coordinate-based ST formats supported by Giotto
- User-provided expression matrices + spatial coordinates + metadata

The pipeline assumes that **data harmonization into Giotto-compatible structures**
happens at ingest time, after which downstream steps are data-source agnostic.

---

## Major Components

### 1) Container Layer

**Location:** `container/`

Responsibilities:
- Provide a fixed R runtime (version-pinned where feasible)
- Install Giotto and core dependencies
- Encapsulate all execution logic for HPC portability

Key files:
- `Singularity.def`: container recipe
- `build.sh`: deterministic image build
- `run.sh`: standardized entrypoint wrapper

Assumptions:
- No outbound internet during runtime
- Read-only container filesystem
- Host directories bind-mounted for inputs/outputs

---

### 2) Orchestration Scripts

**Location:** `scripts/`

Responsibilities:
- Define *what runs* and *in what order*
- Parse configuration
- Call reusable functions
- Write results to disk

Characteristics:
- Executed via `Rscript` inside the container
- No hidden global state
- Minimal logic beyond orchestration

Typical scripts:
- `run_all.R`: end-to-end pipeline
- `qc.R`: quality control only
- `plots.R`: spatial visualization outputs

---

### 3) Reusable R Modules

**Location:** `R/`

Responsibilities:
- Implement pure, reusable functionality
- Encapsulate ST operations independent of execution context

Examples:
- Data validation and filtering
- Normalization and feature selection
- Dimensionality reduction and clustering
- Spatial statistics and plotting helpers

Design principles:
- Explicit inputs and outputs
- Early validation with informative errors
- No filesystem side effects unless explicitly passed paths

---

### 4) Configuration Layer

**Location:** `configs/`

Responsibilities:
- Capture user-configurable parameters
- Decouple pipeline logic from dataset-specific settings

Typical configuration concerns:
- Input data paths
- Output directory
- QC thresholds
- Normalization and clustering parameters
- Plotting options

Configs are treated as **inputs**, not mutable state.

---

### 5) Documentation Layer

**Location:** `docs/`

Responsibilities:
- Preserve project memory and intent

Key documents:
- `decision_log.md`: append-only record of architectural decisions
- `session_notes.md`: append-only daily development notes
- `architecture.md`: this document

---

## Container-as-a-Function Execution Model

The pipeline container is designed to behave as a **pure computational function**.

Conceptually:
outputs = container(inputs, configuration)


### Properties
- A single canonical entrypoint script (e.g., `scripts/run_all.R`)
- Explicit inputs:
  - configuration file
  - data paths
  - runtime options
- Explicit outputs:
  - tables, plots, serialized objects written to a user-defined directory
- No reliance on interactive sessions or container-internal state

### Implications
- The same entrypoint semantics apply whether running locally or in a container.
- Local-first development must preserve function-like behavior.
- The container runtime is an implementation detail, not a user-facing interface.

This model mirrors the execution semantics used in the
`multi-gene-correlations` project and is a core architectural invariant.


---

## Data Flow

User-provided ST data
(expression + spatial coords + metadata)
↓
Ingest & validation
↓
Giotto object
↓
QC & filtering stages
↓
Normalization & feature selection
↓
Dimensionality reduction & clustering
↓
Spatial analysis & visualization
↓
Tables, plots, serialized objects
(written to results/)


At no point is data assumed to be interactive-only or in-memory-only.

---

## Output Conventions

All outputs are written to a predictable directory structure under a user-specified
output path (e.g., `results/`), including:
- QC summaries
- Cluster assignments
- Marker tables
- Spatial plots
- Serialized Giotto objects (RDS)

This ensures:
- reproducibility
- restartability
- downstream reuse

---

## Relationship to Decisions

This architecture reflects decisions logged in `docs/decision_log.md`,
including (but not limited to):
- selection of Giotto as the core ST framework
- Singularity/Apptainer as the container runtime
- script-driven (not notebook-driven) execution

Any future architectural changes should be recorded there before modifying this document.




