FROM rocker/r-ver:4.3.2

ENV DEBIAN_FRONTEND=noninteractive

# ---------------------------------------------------------------------
# Install system dependencies
# ---------------------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        gdebi-core \
        pandoc \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libcairo2-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff5-dev \
        librsvg2-dev \
        libgit2-dev \
        libudunits2-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libsqlite3-dev \
        libglpk-dev \
        libuv1-dev \
        make \
        g++ \
        gfortran && \
    rm -rf /var/lib/apt/lists/*
# ---------------------------------------------------------------------
# Install Quarto (works on both Intel and Apple Silicon)
# ---------------------------------------------------------------------
ARG TARGETARCH

RUN if [ "$TARGETARCH" = "arm64" ]; then \
        QUARTO_ARCH="arm64"; \
    else \
        QUARTO_ARCH="amd64"; \
    fi && \
    wget "https://github.com/quarto-dev/quarto-cli/releases/download/v1.7.34/quarto-1.7.34-linux-${QUARTO_ARCH}.deb" && \
    gdebi --non-interactive "quarto-1.7.34-linux-${QUARTO_ARCH}.deb" && \
    rm "quarto-1.7.34-linux-${QUARTO_ARCH}.deb"

# ---------------------------------------------------------------------
# Install R packages
# ---------------------------------------------------------------------

RUN Rscript -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('readr','dplyr','tidyr','stringr','purrr','tibble','ggplot2','patchwork','scales','knitr','rmarkdown','gt'))"
# ---------------------------------------------------------------------
# Verify installation
# ---------------------------------------------------------------------
RUN quarto check

WORKDIR /work