#!/usr/bin/env bash
set -euo pipefail

OUT=${1:-/tmp/giotto-extra-ca.pem}
shift || true

HOSTS=("$@")
if [[ ${#HOSTS[@]} -eq 0 ]]; then
	HOSTS=(mghp.osn.xsede.org cloud.r-project.org)
fi

tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT

: >"${OUT}"
for host in "${HOSTS[@]}"; do
	echo "Exporting presented TLS chain for ${host}" >&2
	openssl s_client \
		-connect "${host}:443" \
		-servername "${host}" \
		-showcerts </dev/null 2>/dev/null \
		| awk '
			/-----BEGIN CERTIFICATE-----/ { keep = 1 }
			keep { print }
			/-----END CERTIFICATE-----/ { keep = 0 }
		' >"${tmp}"

	if [[ ! -s "${tmp}" ]]; then
		echo "No certificates exported for ${host}" >&2
		exit 1
	fi

	cat "${tmp}" >>"${OUT}"
done

cert_count=$(grep -c 'BEGIN CERTIFICATE' "${OUT}" || true)
echo "Wrote ${cert_count} certificate(s) to ${OUT}" >&2
