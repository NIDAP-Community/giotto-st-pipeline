# Container build and run

This project follows a Container-as-a-Function pattern:

- one stable entrypoint
- explicit bind mounts
- all analysis inputs under `/data`
- all outputs written under `/output`

Pull the published image for normal use, or rebuild it from source on a workstation when you need to modify the container contents.

For this repository, the published GHCR artifact is a lean image:

- `arrow` is excluded during `renv::restore()` by default
- Xenium, h5ad, and Visium inputs using CSV spatial metadata remain supported
- parquet-backed Visium spatial metadata is out of scope for the lean image

This keeps the default build smaller and more reliable for publication and rebuilds.

The source Dockerfile now also provisions a dedicated Linux Python environment for Giotto-backed `h5ad` ingest and clustering helpers. Rebuilds should therefore validate `h5ad` support from inside the image rather than relying on a host Python environment.

## Publish to GHCR

Publish the lean image with one immutable tag plus the floating aliases `lean` and `latest`.
The repository convention is:

- immutable tag: `sha-<git-sha>`
- floating tags: `lean`, `latest`

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

If local Docker rebuilds are blocked by enterprise mirror or signature issues, run the manual GitHub Actions workflow at `.github/workflows/publish-ghcr.yml`. It builds the lean image on GitHub-hosted runners and can publish `sha-<git-sha>`, `lean`, and `latest` without depending on the local workstation network path.

If the workflow fails with a GHCR permission error such as `permission_denied: write_package`, add a repository or organization Actions secret named `GHCR_TOKEN` and rerun the workflow. That token should include at least `write:packages` and `read:packages`; include `repo` as well when the package is tied to a private repository or org policy requires it. The workflow prefers `GHCR_TOKEN` automatically and falls back to `GITHUB_TOKEN` only when no PAT secret is configured.

## Pull the published image

```
docker pull ghcr.io/nidap-community/giotto-st-pipeline:latest
```

For HPC environments, convert the published OCI image directly into a `.sif` artifact:

```
singularity pull giotto-st-pipeline.sif docker://ghcr.io/nidap-community/giotto-st-pipeline:latest
```

## Build the image from source

```
./container/build.sh giotto-st-pipeline:latest
```

The script wraps `docker build` to ensure the context and Dockerfile path are
consistent. Feel free to substitute a different tag as needed.

On Apple Silicon, you can pass extra Docker arguments through the helper:

```
./container/build.sh giotto-st-pipeline:latest --platform linux/amd64
```

If your organization uses an HTTPS inspection proxy or internal root CA, pass it
as a BuildKit secret rather than committing certificates into the repository:

```
./container/build.sh giotto-st-pipeline:latest \
	--secret id=extra_ca,src=/path/to/internal-root-ca.pem
```

The Dockerfile will install that CA into the image before `renv::restore()` runs.

Recommended documentation stance:

- end users should pull the published GHCR image and do not need to manage enterprise CA certificates locally
- maintainers rebuilding from source behind an enterprise TLS-inspecting proxy should inject the internal root CA at build time
- never commit internal CA certificates into this repository

By default, the image excludes the optional `arrow` R package during `renv::restore()` to keep Docker disk usage manageable. This means Visium parquet tissue-position files are not available in the default image build; CSV tissue-position files still work.

On a machine with enough Docker disk space, build the fuller image by clearing the exclude list:

```
./container/build.sh giotto-st-pipeline:latest \
	--build-arg RENV_RESTORE_EXCLUDE= \
	--secret id=extra_ca,src=/path/to/internal-root-ca.pem
```

If this project later publishes both lean and full images, the lean image should remain the default documented rebuild path unless the fuller variant proves equally reliable in CI.

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
