#!/usr/bin/env bash
set -euo pipefail

cd /opt/giotto-st-pipeline

if [[ $# -eq 0 ]]; then
	exec Rscript scripts/run_all.R --help
fi

exec Rscript scripts/run_all.R "$@"
