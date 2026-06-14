#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-ghcr.io/nidap-community/giotto-st-pipeline}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
	cat >&2 <<'USAGE'
Usage: ./container/resolve_image_digest.sh [tag]

Resolve the GHCR registry digest for a published giotto-st-pipeline image tag.
If no tag is supplied, the script uses sha-<current-git-sha> from this checkout.

Examples:
  ./container/resolve_image_digest.sh
  ./container/resolve_image_digest.sh sha-718bf970109a034e01825b3d956a83e3f91ffd9a
  ./container/resolve_image_digest.sh latest

Set IMAGE_NAME to resolve a different image repository.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	usage
	exit 0
fi

TAG="${1:-}"
if [[ -z "$TAG" ]]; then
	TAG="sha-$(git -C "$REPO_ROOT" rev-parse HEAD)"
fi

if [[ "$TAG" == sha256:* ]]; then
	printf '%s\n' "$TAG"
	exit 0
fi

if [[ "$TAG" == *@sha256:* ]]; then
	printf '%s\n' "${TAG##*@}"
	exit 0
fi

if command -v docker >/dev/null 2>&1 && docker buildx version >/dev/null 2>&1; then
	digest="$(docker buildx imagetools inspect "${IMAGE_NAME}:${TAG}" 2>/dev/null | awk '/^Digest:/ { print $2; exit }')"
	if [[ -n "${digest:-}" ]]; then
		printf '%s\n' "$digest"
		exit 0
	fi
fi

if ! command -v curl >/dev/null 2>&1; then
	echo "Could not resolve ${IMAGE_NAME}:${TAG}: install curl or Docker buildx." >&2
	exit 1
fi

repo="${IMAGE_NAME#ghcr.io/}"
token_json="$(curl -fsSL "https://ghcr.io/token?service=ghcr.io&scope=repository:${repo}:pull")"
token="$(printf '%s' "$token_json" | sed -n 's/.*"token":"\([^"]*\)".*/\1/p')"
if [[ -z "$token" ]]; then
	echo "Could not get a GHCR pull token for ${IMAGE_NAME}." >&2
	exit 1
fi

headers="$(curl -fsSLI \
	-H "Authorization: Bearer ${token}" \
	-H "Accept: application/vnd.docker.distribution.manifest.v2+json, application/vnd.oci.image.manifest.v1+json, application/vnd.docker.distribution.manifest.list.v2+json, application/vnd.oci.image.index.v1+json" \
	"https://ghcr.io/v2/${repo}/manifests/${TAG}")"

digest="$(printf '%s\n' "$headers" | awk 'BEGIN { IGNORECASE=1 } /^docker-content-digest:/ { gsub("\r", "", $2); print $2; exit }')"
if [[ -z "$digest" ]]; then
	echo "Could not find Docker-Content-Digest for ${IMAGE_NAME}:${TAG}." >&2
	exit 1
fi

printf '%s\n' "$digest"
