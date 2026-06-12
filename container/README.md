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

Publish the reproducible image with one immutable tag plus floating aliases.
The repository convention is:

- immutable tag: `sha-<full-git-sha>`
- floating tag: `latest` for the full lockfile image

The immutable tag is generated from the committed source revision that contains
the `renv.lock` used by the build. Local publishing refuses to run from a dirty
worktree by default so the source tag cannot accidentally describe a different
set of files than the image contains. It also verifies that the local image's
source revision and `renv.lock` SHA-256 labels match the current checkout before
pushing.

Authenticate Docker to GHCR with a token that includes `write:packages`.
If the package will be private, include `repo` as well.

```
export GHCR_TOKEN=...
echo "$GHCR_TOKEN" | docker login ghcr.io -u <github-username> --password-stdin
./container/publish.sh giotto-st-pipeline:dev
```

Use `--dry-run` to preview tags without pushing:

```
./container/publish.sh giotto-st-pipeline:dev --dry-run
```

Override the immutable tag when you need a release-specific label:

```
./container/publish.sh giotto-st-pipeline:dev --immutable-tag 2026.04.02
```

Only use `--allow-dirty` for explicitly non-reproducible test publishes, and
pair it with a custom immutable tag:

```
./container/publish.sh giotto-st-pipeline:dev \
	--immutable-tag test-2026-06-11 \
	--allow-dirty
```

If local Docker rebuilds are blocked by enterprise mirror or signature issues, run the manual GitHub Actions workflow at `.github/workflows/publish-ghcr.yml`. It builds on GitHub-hosted runners and can publish `sha-<full-git-sha>` and floating tags without depending on the local workstation network path.

If the workflow fails with a GHCR permission error such as `permission_denied: write_package`, add a repository or organization Actions secret named `GHCR_TOKEN` and rerun the workflow. That token should include at least `write:packages` and `read:packages`; include `repo` as well when the package is tied to a private repository or org policy requires it. The workflow prefers `GHCR_TOKEN` automatically and falls back to `GITHUB_TOKEN` only when no PAT secret is configured.

## Pull the published image

```
docker pull ghcr.io/nidap-community/giotto-st-pipeline:latest
```

Use `latest` only for exploratory runs. For a reproducible analysis, pin to the
source commit tag:

```
docker pull ghcr.io/nidap-community/giotto-st-pipeline:sha-<full-git-sha>
```

For byte-identical registry pinning, use the image digest reported by GHCR,
`docker pull`, or `docker inspect`:

```
docker pull ghcr.io/nidap-community/giotto-st-pipeline@sha256:<image-digest>
```

For HPC environments, convert the published OCI image directly into a `.sif` artifact:

```
singularity pull giotto-st-pipeline.sif docker://ghcr.io/nidap-community/giotto-st-pipeline:latest
```

For a pinned HPC artifact, replace `:latest` with `:sha-<full-git-sha>` or
`@sha256:<image-digest>`.

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

## Run the pipeline

Bind (or mount) input and output directories and forward the desired CLI flags
to the `scripts/run_all.R` entrypoint. Example:

```
docker run --rm \
	-v /path/to/xenium:/data:ro \
	-v /path/to/results:/output \
	ghcr.io/nidap-community/giotto-st-pipeline:latest \
	--input_format xenium \
	--input_dir /data/output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834 \
	--output_dir /output/xenium_r1
```

The container entrypoint delegates directly to `scripts/run_all.R`, so all CLI
flags documented for the script are accepted.

## Run with Apptainer/Singularity

Use the helper script to standardize binds in HPC environments:

```
DATA_DIR=$PWD/test_data \
OUTPUT_DIR=$PWD/test_output \
./container/run_apptainer.sh /path/to/giotto-st-pipeline.sif \
	--input_format h5ad \
	--input_path /data/sample.h5ad \
	--output_dir /output/sample_h5ad \
	--project_id sample_h5ad
```

`DATA_DIR`, `OUTPUT_DIR`, and optional `EXTRA_BINDS` control the bind mounts.
The helper also exports `R_LIBS_USER=/usr/local/lib/R/site-library` so packaged
R dependencies remain discoverable when `HOME` is remapped by the runtime.

## Exporting for Apptainer/Singularity

After building a local Docker image you can export it to a tarball and convert it into a `.sif` artifact:

```
docker save giotto-st-pipeline:latest -o giotto-st-pipeline.tar
singularity build giotto-st-pipeline.sif docker-archive://giotto-st-pipeline.tar
```

`giotto-st-pipeline.sif` can then be copied to the HPC environment and executed
with `apptainer run` or `singularity run`, binding the same input/output
directories as shown above.
