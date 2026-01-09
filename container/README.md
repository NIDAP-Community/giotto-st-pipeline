# Container build and run

This project uses Docker as the reference build target so the resulting image can
be exported and executed with Apptainer/Singularity on clusters that do not
support Docker directly. All examples below assume you are in the repository
root on a workstation that has Docker available.

## Build the image

```
./container/build.sh giotto-st-pipeline:latest
```

The script wraps `docker build` to ensure the context and Dockerfile path are
consistent. Feel free to substitute a different tag as needed.

## Run the pipeline

Bind (or mount) input and output directories and forward the desired CLI flags
to the `scripts/run_all.R` entrypoint. Example:

```
docker run --rm \
	-v /path/to/xenium:/input:ro \
	-v /path/to/results:/output \
	giotto-st-pipeline:latest \
	--input_format xenium \
	--input_dir /input/output-XETG00202__0024834_Right__SCAF04264_Right_R1__20240912__162834 \
	--output_dir /output/xenium_r1
```

The container entrypoint delegates directly to `scripts/run_all.R`, so all CLI
flags documented for the script are accepted.

## Exporting for Apptainer/Singularity

After building the Docker image you can push it to a registry or export it to a
local tarball and convert it into a `.sif` artifact:

```
docker save giotto-st-pipeline:latest -o giotto-st-pipeline.tar
singularity build giotto-st-pipeline.sif docker-archive://giotto-st-pipeline.tar
```

`giotto-st-pipeline.sif` can then be copied to the HPC environment and executed
with `apptainer run` or `singularity run`, binding the same input/output
directories as shown above.
