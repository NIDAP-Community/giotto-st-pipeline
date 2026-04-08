# syntax=docker/dockerfile:1.6

FROM rocker/r-ver:4.4.3

LABEL org.opencontainers.image.title="giotto-st-pipeline" \
    org.opencontainers.image.description="Giotto spatial transcriptomics workflow with a staged Container-as-a-Function entrypoint" \
    org.opencontainers.image.url="https://github.com/NIDAP-Community/giotto-st-pipeline" \
    org.opencontainers.image.source="https://github.com/NIDAP-Community/giotto-st-pipeline" \
    org.opencontainers.image.documentation="https://github.com/NIDAP-Community/giotto-st-pipeline#readme"

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_CACHE_ENABLED=TRUE \
    RETICULATE_PYTHON=/opt/giotto-python/bin/python

ARG RENV_RESTORE_EXCLUDE=arrow

WORKDIR /opt/giotto-st-pipeline

RUN sed -i 's|http://ports.ubuntu.com/ubuntu-ports|https://ports.ubuntu.com/ubuntu-ports|g; s|http://archive.ubuntu.com/ubuntu|https://archive.ubuntu.com/ubuntu|g' /etc/apt/sources.list /etc/apt/sources.list.d/*.sources 2>/dev/null || true

RUN --mount=type=secret,id=extra_ca,target=/run/secrets/extra_ca,required=false \
    if [ -f /run/secrets/extra_ca ]; then \
    cp /run/secrets/extra_ca /usr/local/share/ca-certificates/extra-ca.crt && \
    update-ca-certificates; \
    fi

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    g++ \
    gfortran \
    git \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgdal-dev \
    libgeos-dev \
    libglpk-dev \
    libgsl-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libjpeg-dev \
    libmagick++-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libudunits2-dev \
    libxml2-dev \
    pkg-config \
    python3 \
    python3-venv \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/giotto-python \
    && /opt/giotto-python/bin/pip install --no-cache-dir --upgrade pip setuptools wheel \
    && /opt/giotto-python/bin/pip install --no-cache-dir \
    "anndata>=0.7" \
    scipy \
    python-igraph \
    leidenalg \
    networkx \
    scikit-learn \
    python-louvain

# Copy renv metadata first to maximise Docker layer caching
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e 'install.packages("renv", repos = "https://cloud.r-project.org"); exclude <- Sys.getenv("RENV_RESTORE_EXCLUDE", unset = ""); if (nzchar(exclude)) { renv::restore(confirm = FALSE, exclude = strsplit(exclude, ",", fixed = TRUE)[[1]]) } else { renv::restore(confirm = FALSE) }'

COPY .Rprofile .Rprofile

# Copy the remainder of the repository (scripts, configs, docs)
COPY . .

RUN chmod +x container/run.sh container/build.sh

ENTRYPOINT ["/opt/giotto-st-pipeline/container/run.sh"]
CMD ["--help"]
