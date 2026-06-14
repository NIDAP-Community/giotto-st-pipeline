# Container build and run

This project follows a Container-as-a-Function pattern:

- one stable entrypoint
- explicit bind mounts
- all analysis inputs under `/data`
- all outputs written under `/output`

Pull the published image for normal use, or rebuild it from source on a workstation when you need to modify the container contents.

For this repository, the reproducible GHCR artifact should restore the full
locked runtime:

- the Docker base image is pinned by digest
- R packages are restored from `renv.lock`
- the `renv` bootstrap version is pinned before restore
- Python packages are installed from pinned requirement files
- Xenium, h5ad, and Visium inputs remain supported

This makes the default build heavier, but it avoids silently publishing an image
that omits packages present in the lockfile. Build a lean variant only when you
explicitly accept the reduced feature surface.

The source Dockerfile now also provisions a dedicated Linux Python environment for Giotto-backed `h5ad` ingest and clustering helpers. Rebuilds should therefore validate `h5ad` support from inside the image rather than relying on a host Python environment.

## Publish to GHCR

Publication uses the GitHub Actions workflow at `.github/workflows/publish-ghcr.yml`.
The workflow reads `VERSION` and publishes three tags per release:

- immutable tag: `sha-<full-git-sha>`
- version tag: `<version>` (e.g., `0.1.0`)
- floating tag: `latest`

To publish, trigger the workflow manually from the Actions tab (or via `gh`):

```bash
gh workflow run publish-ghcr.yml
```

If the workflow fails with a GHCR permission error such as `permission_denied: write_package`, add a repository or organization Actions secret named `GHCR_TOKEN` and rerun the workflow. That token should include at least `write:packages` and `read:packages`; include `repo` as well when the package is tied to a private repository or org policy requires it. The workflow prefers `GHCR_TOKEN` automatically and falls back to `GITHUB_TOKEN` only when no PAT secret is configured.

## Run With Apptainer/Singularity

Most HPC users should pull the published OCI image into a `.sif` artifact and
run that file on a compute node. Pull by version tag for reproducibility:

```
module load singularity

IMAGE_REPO="ghcr.io/nidap-community/giotto-st-pipeline"
IMAGE_TAG="0.1.0"

singularity pull giotto-st-pipeline.sif "docker://${IMAGE_REPO}:${IMAGE_TAG}"
```

For bit-level pinning, resolve the digest with `skopeo` (if available):

```
module load skopeo

IMAGE_DIGEST="$(skopeo inspect --format '{{.Digest}}' "docker://${IMAGE_REPO}:${IMAGE_TAG}")"
IMAGE_URI="${IMAGE_REPO}@${IMAGE_DIGEST}"

singularity pull giotto-st-pipeline.sif "docker://${IMAGE_URI}"
```

Bind input data at `/data`, output at `/output`, and pass a config file whose
paths match those container mounts:

```
singularity run --cleanenv \
	--bind /path/to/spaceranger/outs:/data:ro \
	--bind /path/to/results:/output \
	--bind "$PWD/configs":/configs:ro \
	giotto-st-pipeline.sif \
	--config /configs/my_visium_container.yaml
```

Use the helper script to standardize binds in HPC environments and pass SIF
provenance into the run metadata:

```
DATA_DIR=/path/to/spaceranger/outs \
OUTPUT_DIR=/path/to/results \
EXTRA_BINDS="$PWD/configs:/configs:ro" \
CONTAINER_IMAGE_REF="$IMAGE_URI" \
CONTAINER_IMAGE_DIGEST="${IMAGE_DIGEST:-}" \
./container/run_apptainer.sh giotto-st-pipeline.sif \
	--config /configs/my_visium_container.yaml
```

`DATA_DIR`, `OUTPUT_DIR`, and optional `EXTRA_BINDS` control the bind mounts.
The helper exports the runtime name, SIF path, SIF SHA-256, and optional source
image reference so `metadata/provenance.json` captures what was executed.

## Pull the published image with Docker

Docker is useful for local workstation runs and image maintenance. Pull by
version tag:

```
docker pull ghcr.io/nidap-community/giotto-st-pipeline:0.1.0
```

Inspect the version and source commit baked into the image:

```
IMAGE_REPO="ghcr.io/nidap-community/giotto-st-pipeline"
IMAGE_TAG="0.1.0"
docker inspect "${IMAGE_REPO}:${IMAGE_TAG}" \
	--format 'version={{ index .Config.Labels "org.opencontainers.image.version" }} commit={{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

## Build the image from source

```
./container/build.sh giotto-st-pipeline:latest
```

The script wraps `docker build` to ensure the context and Dockerfile path are
consistent. Feel free to substitute a different tag as needed.

Local builds stamp the image with `org.opencontainers.image.revision` and
`org.opencontainers.image.renv-lock-sha256` labels. Inspect them with:

```
docker inspect giotto-st-pipeline:latest \
	--format '{{ index .Config.Labels "org.opencontainers.image.revision" }} {{ index .Config.Labels "org.opencontainers.image.renv-lock-sha256" }}'
```

On Apple Silicon, you can pass extra Docker arguments through the helper:

```
./container/build.sh giotto-st-pipeline:latest --platform linux/amd64
```

The default build keeps the Ubuntu apt sources from the upstream Rocker image.
If your organization requires HTTPS apt sources and uses an HTTPS inspection
proxy or internal root CA, enable the HTTPS rewrite and pass the CA as a BuildKit
secret rather than committing certificates into the repository:

```
./container/build.sh giotto-st-pipeline:latest \
	--build-arg APT_USE_HTTPS=true \
	--build-arg EXTRA_CA_CACHE_BUST="$(date +%s)" \
	--secret id=extra_ca,src=/path/to/internal-root-ca.pem
```

The Dockerfile will install that CA into the image before `renv::restore()` runs.
The `EXTRA_CA_CACHE_BUST` argument forces Docker to refresh the CA-install layer
when the secret changes; BuildKit does not include secret contents in the cache
key.

If Docker fails with `SSL certificate problem: self-signed certificate in
certificate chain` while querying Bioconductor, the failing host may be a
redirect target rather than `bioconductor.org` itself. In NIH network
environments, Bioconductor 3.20 package indexes can redirect to
`mghp.osn.xsede.org`, which may be presented through the NIH inspection CA.
Generate a temporary PEM from the currently presented TLS chains, then pass it
as the same BuildKit secret:

```
./container/export_tls_chain_ca.sh /tmp/giotto-extra-ca.pem

./container/build.sh giotto-st-pipeline:latest \
	--build-arg EXTRA_CA_CACHE_BUST="$(date +%s)" \
	--secret id=extra_ca,src=/tmp/giotto-extra-ca.pem
```

Treat `/tmp/giotto-extra-ca.pem` as a local build artifact. Do not commit it.

Recommended documentation stance:

- end users should pull the published GHCR image and do not need to manage enterprise CA certificates locally
- maintainers rebuilding from source behind an enterprise TLS-inspecting proxy should inject the internal root CA at build time
- never commit internal CA certificates into this repository

By default, the image restores the full `renv.lock`, including optional packages
such as `arrow`. If you need a smaller local debugging image, make that tradeoff
explicit:

```
./container/build.sh giotto-st-pipeline:lean \
	--build-arg RENV_RESTORE_EXCLUDE=arrow
```

This lean variant keeps Xenium, h5ad, and Visium inputs using CSV spatial
metadata, but parquet-backed Visium spatial metadata is out of scope.

The default build uses `http://cloud.r-project.org` as the CRAN mirror so local Docker Desktop builds work in environments where HTTPS package mirrors are intercepted by an enterprise proxy. Override it when your build environment has normal HTTPS CRAN access:

```
./container/build.sh giotto-st-pipeline:latest \
	--build-arg CRAN_REPO=https://cloud.r-project.org
```

The base image is pinned by digest. To refresh it intentionally, update the
`FROM rocker/r-ver:4.4.3@sha256:...` line in `container/Dockerfile` and rebuild
from a clean checkout.

R packages are controlled by `renv.lock`. Python packages are controlled by:

- `container/python-bootstrap-requirements.txt`
- `container/python-requirements.txt`

When changing dependencies, update the appropriate lockfile or requirement file
in the same commit as the code that needs the change.

Ubuntu apt packages are still resolved from the upstream Ubuntu repositories at
build time. For strict archival reproducibility, build against a pinned Ubuntu
snapshot or an internal frozen mirror and record that snapshot with the image
tag. The digest-pinned base image controls the starting filesystem, but it does
not freeze packages installed after `apt-get update`.

To rebuild the full image while also injecting an enterprise CA:

```
./container/build.sh giotto-st-pipeline:latest \
	--secret id=extra_ca,src=/path/to/internal-root-ca.pem
```

When local rebuilds fail because package mirrors or CA chains behave differently on the workstation, prefer the GitHub Actions workflow for publication and use local Docker builds only for iterative debugging.

## Bind mount layout

The container expects explicit mounts:

- `/data`: read-only inputs such as Xenium directories, Visium outputs, `.h5ad` files, and configs
- `/output`: writable results directory

Always mount host paths into those locations when running the image.

## Run the pipeline with Docker

Bind (or mount) input and output directories and forward the desired CLI flags
to the `scripts/run_all.R` entrypoint. Example:

```
docker run --rm \
	-v /path/to/xenium:/data:ro \
	-v /path/to/results:/output \
	-v "$PWD/configs":/configs:ro \
	ghcr.io/nidap-community/giotto-st-pipeline:latest \
	--config /configs/my_xenium_container.yaml
```

The container entrypoint delegates directly to `scripts/run_all.R`, so all CLI
flags documented for the script are accepted.

## Exporting for Apptainer/Singularity

After building a local Docker image you can export it to a tarball and convert it into a `.sif` artifact:

```
docker save giotto-st-pipeline:latest -o giotto-st-pipeline.tar
singularity build giotto-st-pipeline.sif docker-archive://giotto-st-pipeline.tar
```

`giotto-st-pipeline.sif` can then be copied to the HPC environment and executed
with `apptainer run` or `singularity run`, binding the same input/output
directories as shown above.
