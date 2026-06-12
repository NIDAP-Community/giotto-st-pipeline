#!/usr/bin/env bash
set -euo pipefail

LOCAL_IMAGE=${1:-giotto-st-pipeline:dev}
if [[ $# -gt 0 ]]; then
	shift
fi

IMMUTABLE_TAG=""
PUBLISH_FLOATING_TAGS=true
DRY_RUN=false
ALLOW_DIRTY=false

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
		--allow-dirty)
			ALLOW_DIRTY=true
			shift
			;;
		*)
			echo "Usage: ./container/publish.sh [local-image] [--immutable-tag <tag>] [--no-floating-tags] [--dry-run] [--allow-dirty]" >&2
			exit 1
			;;
	esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

GHCR_OWNER=${GHCR_OWNER:-nidap-community}
GHCR_IMAGE_NAME=${GHCR_IMAGE_NAME:-giotto-st-pipeline}
REMOTE_BASE="ghcr.io/${GHCR_OWNER}/${GHCR_IMAGE_NAME}"
CURRENT_GIT_SHA=$(git -C "${REPO_ROOT}" rev-parse HEAD)
CURRENT_RENV_LOCK_SHA256=$(shasum -a 256 "${REPO_ROOT}/renv.lock" | awk '{ print $1 }')

if [[ -z "${IMMUTABLE_TAG}" ]]; then
	IMMUTABLE_TAG="sha-${CURRENT_GIT_SHA}"
fi

if [[ "${ALLOW_DIRTY}" != true ]] && [[ -n "$(git -C "${REPO_ROOT}" status --porcelain)" ]]; then
	echo "Refusing to publish from a dirty worktree." >&2
	echo "Commit the exact source state first so ${IMMUTABLE_TAG} matches the source and renv.lock used to build the image." >&2
	echo "For an explicitly non-reproducible test publish, rerun with --allow-dirty and a custom --immutable-tag." >&2
	exit 1
fi

LOCAL_IMAGE_ID=$(docker images --no-trunc --format '{{.Repository}}:{{.Tag}} {{.ID}}' | awk -v ref="${LOCAL_IMAGE}" '$1 == ref { print $2; exit }')

if [[ -z "${LOCAL_IMAGE_ID}" ]]; then
	echo "Local image not found: ${LOCAL_IMAGE}" >&2
	echo "Build it first with ./container/build.sh ${LOCAL_IMAGE}" >&2
	exit 1
fi

if [[ "${ALLOW_DIRTY}" != true ]]; then
	IMAGE_SOURCE_COMMIT=$(docker image inspect --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}' "${LOCAL_IMAGE_ID}" 2>/dev/null || true)
	IMAGE_RENV_LOCK_SHA256=$(docker image inspect --format '{{ index .Config.Labels "org.opencontainers.image.renv-lock-sha256" }}' "${LOCAL_IMAGE_ID}" 2>/dev/null || true)
	if [[ "${IMAGE_SOURCE_COMMIT}" != "${CURRENT_GIT_SHA}" ]] || [[ "${IMAGE_RENV_LOCK_SHA256}" != "${CURRENT_RENV_LOCK_SHA256}" ]]; then
		echo "Refusing to publish: local image metadata does not match the current committed source." >&2
		echo "Expected source revision: ${CURRENT_GIT_SHA}" >&2
		echo "Image source revision:    ${IMAGE_SOURCE_COMMIT:-<missing>}" >&2
		echo "Expected renv.lock SHA:   ${CURRENT_RENV_LOCK_SHA256}" >&2
		echo "Image renv.lock SHA:      ${IMAGE_RENV_LOCK_SHA256:-<missing>}" >&2
		echo "Rebuild with ./container/build.sh ${LOCAL_IMAGE}, then publish again." >&2
		exit 1
	fi
fi

tags=("${IMMUTABLE_TAG}")
if [[ "${PUBLISH_FLOATING_TAGS}" == true ]]; then
	tags+=("latest")
fi

echo "Local image: ${LOCAL_IMAGE}"
echo "Resolved image ID: ${LOCAL_IMAGE_ID}"
echo "Remote image: ${REMOTE_BASE}"
printf 'Tags to publish:\n'
for tag in "${tags[@]}"; do
	echo "  - ${tag}"
done

if [[ "${DRY_RUN}" == true ]]; then
	exit 0
fi

for tag in "${tags[@]}"; do
	docker tag "${LOCAL_IMAGE_ID}" "${REMOTE_BASE}:${tag}"
done

for tag in "${tags[@]}"; do
	docker push "${REMOTE_BASE}:${tag}"
done
