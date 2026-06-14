#!/usr/bin/env bash
set -euo pipefail

IMAGE_TAG=${1:-giotto-st-pipeline:latest}
if [[ $# -gt 0 ]]; then
	shift
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SOURCE_COMMIT="$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
RENV_LOCK_SHA256="$(shasum -a 256 "${REPO_ROOT}/renv.lock" | awk '{ print $1 }')"
PIPELINE_VERSION="$(cat "${REPO_ROOT}/VERSION" 2>/dev/null || echo unknown)"

docker build \
	"$@" \
	--build-arg "SOURCE_COMMIT=${SOURCE_COMMIT}" \
	--build-arg "RENV_LOCK_SHA256=${RENV_LOCK_SHA256}" \
	--build-arg "PIPELINE_VERSION=${PIPELINE_VERSION}" \
	--tag "${IMAGE_TAG}" \
	--file "${SCRIPT_DIR}/Dockerfile" \
	"${REPO_ROOT}"
