# syntax=docker/dockerfile:1.6

FROM rocker/r-ver:4.4.3

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_CACHE_ENABLED=TRUE

WORKDIR /opt/giotto-st-pipeline

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
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
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy renv metadata first to maximise Docker layer caching
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY renv/settings.dcf renv/settings.dcf
COPY .Rprofile .Rprofile

RUN R -e 'install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore(confirm = FALSE)'

# Copy the remainder of the repository (scripts, configs, docs)
COPY . .

RUN chmod +x container/run.sh container/build.sh

ENTRYPOINT ["/opt/giotto-st-pipeline/container/run.sh"]
CMD ["--help"]
