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
CONTAINER_IMAGE_REF="${CONTAINER_IMAGE_REF:-}"
CONTAINER_IMAGE_DIGEST="${CONTAINER_IMAGE_DIGEST:-}"

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

if command -v sha256sum >/dev/null 2>&1; then
	SIF_SHA256="$(sha256sum "$IMAGE_PATH" | awk '{ print $1 }')"
elif command -v shasum >/dev/null 2>&1; then
	SIF_SHA256="$(shasum -a 256 "$IMAGE_PATH" | awk '{ print $1 }')"
else
	SIF_SHA256=""
fi

export SINGULARITYENV_GIOTTO_PIPELINE_CONTAINER_RUNTIME="$RUNTIME"
export APPTAINERENV_GIOTTO_PIPELINE_CONTAINER_RUNTIME="$RUNTIME"
export SINGULARITYENV_GIOTTO_PIPELINE_SIF_PATH="$IMAGE_PATH"
export APPTAINERENV_GIOTTO_PIPELINE_SIF_PATH="$IMAGE_PATH"
export SINGULARITYENV_GIOTTO_PIPELINE_SIF_SHA256="$SIF_SHA256"
export APPTAINERENV_GIOTTO_PIPELINE_SIF_SHA256="$SIF_SHA256"

if [[ -n "$CONTAINER_IMAGE_REF" ]]; then
	export SINGULARITYENV_GIOTTO_PIPELINE_CONTAINER_IMAGE="$CONTAINER_IMAGE_REF"
	export APPTAINERENV_GIOTTO_PIPELINE_CONTAINER_IMAGE="$CONTAINER_IMAGE_REF"
fi

if [[ -n "$CONTAINER_IMAGE_DIGEST" ]]; then
	export SINGULARITYENV_GIOTTO_PIPELINE_CONTAINER_DIGEST="$CONTAINER_IMAGE_DIGEST"
	export APPTAINERENV_GIOTTO_PIPELINE_CONTAINER_DIGEST="$CONTAINER_IMAGE_DIGEST"
fi

R_LIBS_USER="${R_LIBS_USER:-/usr/local/lib/R/site-library}" \
	exec "$RUNTIME" run "${bind_args[@]}" "$IMAGE_PATH" "$@"
