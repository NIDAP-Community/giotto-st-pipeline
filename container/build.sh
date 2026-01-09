#!/usr/bin/env bash
set -euo pipefail

IMAGE_TAG=${1:-giotto-st-pipeline:latest}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

docker build \
	--tag "${IMAGE_TAG}" \
	--file "${SCRIPT_DIR}/Dockerfile" \
	"${REPO_ROOT}"
