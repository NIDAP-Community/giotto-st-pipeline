#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
	echo "Usage: DATA_DIR=/host/data OUTPUT_DIR=/host/output $0 <image.sif> [pipeline args...]" >&2
	exit 1
fi

if command -v apptainer >/dev/null 2>&1; then
	RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
	RUNTIME="singularity"
else
	echo "Neither apptainer nor singularity is available on PATH." >&2
	exit 1
fi

IMAGE_PATH="$1"
shift

DATA_DIR="${DATA_DIR:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
EXTRA_BINDS="${EXTRA_BINDS:-}"

if [[ -z "$DATA_DIR" ]]; then
	echo "Set DATA_DIR to the host directory that should be mounted at /data." >&2
	exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
	echo "Set OUTPUT_DIR to the host directory that should be mounted at /output." >&2
	exit 1
fi

mkdir -p "$OUTPUT_DIR"

bind_args=(
	"--bind" "$DATA_DIR:/data"
	"--bind" "$OUTPUT_DIR:/output"
)

if [[ -n "$EXTRA_BINDS" ]]; then
	IFS=',' read -r -a extra_bind_list <<< "$EXTRA_BINDS"
	for bind_spec in "${extra_bind_list[@]}"; do
		[[ -n "$bind_spec" ]] || continue
		bind_args+=("--bind" "$bind_spec")
	done
fi

R_LIBS_USER="${R_LIBS_USER:-/usr/local/lib/R/site-library}" \
	exec "$RUNTIME" run "${bind_args[@]}" "$IMAGE_PATH" "$@"
