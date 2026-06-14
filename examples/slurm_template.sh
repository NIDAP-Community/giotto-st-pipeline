#!/bin/bash
#SBATCH --job-name=giotto-st
#SBATCH --output=giotto-st_%j.out
#SBATCH --error=giotto-st_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=norm

# ===========================================================================
# Giotto ST Pipeline — SLURM Job Template
# ===========================================================================
# Usage:
#   1. Copy this file and edit the variables below.
#   2. Submit with: sbatch my_job.sh
#
# Prerequisites:
#   - A .sif container image (see QUICKSTART.md §3)
#   - A YAML config file with container-facing paths (/data, /output)
# ===========================================================================

# ---- User-editable variables ------------------------------------------------

# Path to the .sif file on shared storage
SIF_IMAGE="/path/to/giotto-st-pipeline.sif"

# Host directory containing your ST data (mounted read-only at /data)
DATA_DIR="/path/to/spaceranger/outs"

# Host directory for results (mounted read-write at /output)
OUTPUT_DIR="/path/to/results"

# Config file path on host (will be mounted at /configs inside container)
CONFIG_DIR="$PWD/configs"
CONFIG_FILE="/configs/my_visium_container.yaml"

# Optional: extra pipeline flags (e.g., --cluster_resolution 0.8)
EXTRA_ARGS=""

# ---- End user-editable section ----------------------------------------------

set -euo pipefail

module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true

if ! command -v singularity &>/dev/null && ! command -v apptainer &>/dev/null; then
    echo "ERROR: Neither singularity nor apptainer found. Load the appropriate module." >&2
    exit 1
fi

RUNTIME="singularity"
command -v apptainer &>/dev/null && RUNTIME="apptainer"

mkdir -p "$OUTPUT_DIR"

echo "=============================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "Started:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "SIF image:    $SIF_IMAGE"
echo "Data dir:     $DATA_DIR"
echo "Output dir:   $OUTPUT_DIR"
echo "Config:       $CONFIG_FILE"
echo "Runtime:      $RUNTIME"
echo "=============================================="

# Record container provenance
if command -v sha256sum &>/dev/null; then
    sha256sum "$SIF_IMAGE" > "${OUTPUT_DIR}/giotto-st-pipeline.sif.sha256"
elif command -v shasum &>/dev/null; then
    shasum -a 256 "$SIF_IMAGE" > "${OUTPUT_DIR}/giotto-st-pipeline.sif.sha256"
fi

"$RUNTIME" run --cleanenv \
    --bind "${DATA_DIR}:/data:ro" \
    --bind "${OUTPUT_DIR}:/output" \
    --bind "${CONFIG_DIR}:/configs:ro" \
    "$SIF_IMAGE" \
    --config "$CONFIG_FILE" \
    $EXTRA_ARGS

echo "Completed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
