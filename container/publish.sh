#!/usr/bin/env bash
set -euo pipefail

LOCAL_IMAGE=${1:-giotto-st-pipeline:dev}
if [[ $# -gt 0 ]]; then
	shift
fi

IMMUTABLE_TAG=""
PUBLISH_FLOATING_TAGS=true
DRY_RUN=false

while [[ $# -gt 0 ]]; do
	case "$1" in
		--immutable-tag)
			IMMUTABLE_TAG=${2:?missing value for --immutable-tag}
			shift 2
			;;
		--no-floating-tags)
			PUBLISH_FLOATING_TAGS=false
			shift
			;;
		--dry-run)
			DRY_RUN=true
			shift
			;;
		*)
			echo "Usage: ./container/publish.sh [local-image] [--immutable-tag <tag>] [--no-floating-tags] [--dry-run]" >&2
			exit 1
			;;
	esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

GHCR_OWNER=${GHCR_OWNER:-nidap-community}
GHCR_IMAGE_NAME=${GHCR_IMAGE_NAME:-giotto-st-pipeline}
REMOTE_BASE="ghcr.io/${GHCR_OWNER}/${GHCR_IMAGE_NAME}"

if [[ -z "${IMMUTABLE_TAG}" ]]; then
	GIT_SHA=$(git -C "${REPO_ROOT}" rev-parse --short HEAD)
	IMMUTABLE_TAG="sha-${GIT_SHA}"
fi

LOCAL_IMAGE_ID=$(docker images --no-trunc --format '{{.Repository}}:{{.Tag}} {{.ID}}' | awk -v ref="${LOCAL_IMAGE}" '$1 == ref { print $2; exit }')

if [[ -z "${LOCAL_IMAGE_ID}" ]]; then
	echo "Local image not found: ${LOCAL_IMAGE}" >&2
	echo "Build it first with ./container/build.sh ${LOCAL_IMAGE}" >&2
	exit 1
fi

tags=("${IMMUTABLE_TAG}")
if [[ "${PUBLISH_FLOATING_TAGS}" == true ]]; then
	tags+=("lean" "latest")
fi

echo "Local image: ${LOCAL_IMAGE}"
echo "Resolved image ID: ${LOCAL_IMAGE_ID}"
echo "Remote image: ${REMOTE_BASE}"
printf 'Tags to publish:\n'
for tag in "${tags[@]}"; do
	echo "  - ${tag}"
	docker tag "${LOCAL_IMAGE_ID}" "${REMOTE_BASE}:${tag}"
done

if [[ "${DRY_RUN}" == true ]]; then
	exit 0
fi

for tag in "${tags[@]}"; do
	docker push "${REMOTE_BASE}:${tag}"
done