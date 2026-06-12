# syntax=docker/dockerfile:1.6

FROM rocker/r-ver:4.4.3@sha256:3dae5d2eeddf74f10e0a81fb6b7ae350295e288000304f438b844b2c1e00fe2c

ARG SOURCE_COMMIT=unknown
ARG RENV_LOCK_SHA256=unknown

LABEL org.opencontainers.image.title="giotto-st-pipeline" \
    org.opencontainers.image.description="Giotto spatial transcriptomics workflow with a staged Container-as-a-Function entrypoint" \
    org.opencontainers.image.url="https://github.com/NIDAP-Community/giotto-st-pipeline" \
    org.opencontainers.image.source="https://github.com/NIDAP-Community/giotto-st-pipeline" \
    org.opencontainers.image.documentation="https://github.com/NIDAP-Community/giotto-st-pipeline#readme" \
    org.opencontainers.image.revision="${SOURCE_COMMIT}" \
    org.opencontainers.image.renv-lock-sha256="${RENV_LOCK_SHA256}"

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_CACHE_ENABLED=TRUE \
    RETICULATE_PYTHON=/opt/giotto-python/bin/python \
    SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt \
    CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

ARG RENV_RESTORE_EXCLUDE=
ARG RENV_VERSION=1.1.5
ARG APT_USE_HTTPS=false
ARG CRAN_REPO=http://cloud.r-project.org
ARG EXTRA_CA_CACHE_BUST=0

ENV RENV_CONFIG_REPOS_OVERRIDE=${CRAN_REPO}

WORKDIR /opt/giotto-st-pipeline

RUN if [ "${APT_USE_HTTPS}" = "true" ]; then \
    sed -i 's|http://ports.ubuntu.com/ubuntu-ports|https://ports.ubuntu.com/ubuntu-ports|g; s|http://archive.ubuntu.com/ubuntu|https://archive.ubuntu.com/ubuntu|g' /etc/apt/sources.list /etc/apt/sources.list.d/*.sources 2>/dev/null || true; \
    fi

RUN --mount=type=secret,id=extra_ca,target=/run/secrets/extra_ca,required=false \
    if [ -f /run/secrets/extra_ca ]; then \
    cp /run/secrets/extra_ca /usr/local/share/ca-certificates/extra-ca.crt && \
    update-ca-certificates && \
    cat /run/secrets/extra_ca >> /etc/ssl/certs/ca-certificates.crt; \
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
    libgmp3-dev \
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
    pandoc \
    pkg-config \
    python3 \
    python3-venv \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

COPY container/python-bootstrap-requirements.txt container/python-bootstrap-requirements.txt
COPY container/python-requirements.txt container/python-requirements.txt

RUN python3 -m venv /opt/giotto-python \
    && /opt/giotto-python/bin/pip install --no-cache-dir --requirement container/python-bootstrap-requirements.txt \
    && /opt/giotto-python/bin/pip install --no-cache-dir --no-deps --requirement container/python-requirements.txt \
    && /opt/giotto-python/bin/pip check

RUN --mount=type=secret,id=extra_ca,target=/run/secrets/extra_ca,required=false \
    echo "extra CA cache bust: ${EXTRA_CA_CACHE_BUST}" >/dev/null && \
    if [ -f /run/secrets/extra_ca ]; then \
    cp /run/secrets/extra_ca /usr/local/share/ca-certificates/extra-ca.crt && \
    update-ca-certificates && \
    cat /run/secrets/extra_ca >> /etc/ssl/certs/ca-certificates.crt; \
    fi

# Copy renv metadata first to maximise Docker layer caching
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e 'cran <- sub("/+$", "", Sys.getenv("CRAN_REPO", unset = Sys.getenv("RENV_CONFIG_REPOS_OVERRIDE", unset = "http://cloud.r-project.org"))); renv_version <- Sys.getenv("RENV_VERSION", unset = "1.1.5"); options(repos = c(CRAN = cran)); install.packages(sprintf("%s/src/contrib/Archive/renv/renv_%s.tar.gz", cran, renv_version), repos = NULL, type = "source"); exclude <- Sys.getenv("RENV_RESTORE_EXCLUDE", unset = ""); if (nzchar(exclude)) { renv::restore(confirm = FALSE, exclude = strsplit(exclude, ",", fixed = TRUE)[[1]]) } else { renv::restore(confirm = FALSE) }; renv::status()'

COPY .Rprofile .Rprofile

# Copy the remainder of the repository (scripts, configs, docs)
COPY . .

RUN chmod +x container/run.sh container/build.sh

ENTRYPOINT ["/opt/giotto-st-pipeline/container/run.sh"]
CMD ["--help"]
