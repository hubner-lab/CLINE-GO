# Base image with R 4.5 (requires Bioconductor 3.22)
FROM rocker/shiny:4.5

# Permissions
RUN mkdir -p /.cache && chmod -R 777 /.cache

# Metadata
LABEL maintainer="potapgene@gmail.com"

# Update the system and install CLI tools
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    pipx \
    bash \
    wget \
    curl \
    git \
    gfortran \
    build-essential \
    zlib1g-dev \
    pkg-config \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libglpk40 \
    libuv1 \
    bcftools \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install plink1.9
ENV PLINK_VERSION=20231211
RUN mkdir -p /tmp/plink && \
  cd /tmp/plink && \
  curl -fsSL -o plink.zip "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK_VERSION}.zip" && \
  unzip plink.zip && \
  mv plink /bin/plink && \
  cd $HOME && \
  rm -rf /tmp/plink

# Configure pipx and Install Snakemake
ENV PIPX_HOME=/usr/local/pipx
ENV PIPX_BIN_DIR=/usr/local/bin
RUN pipx install snakemake && pipx inject snakemake pulp==2.7

# Python packages for WZA block-level p-value aggregation
RUN pip3 install --break-system-packages numpy==1.26.4 pandas==2.2.2 scipy==1.13.1

# Set working directory
WORKDIR /pipeline

# Install vcftools
RUN set -ex \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar zxf vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
    && ./configure --prefix=/app/vcftools-0.1.16 \
    && make \
    && make install \
    && cd .. && rm -rf vcftools-0.1.16*

#=============================================================================
# R PACKAGES WITH PINNED VERSIONS
# Last verified: 2025-01-25
# Bioconductor: 3.22 (for R 4.5)
#=============================================================================

# Install remotes first for version-pinned installations
RUN Rscript -e "install.packages('remotes')"

# Core tidyverse and data manipulation (pinned versions)
RUN Rscript -e " \
    remotes::install_version('ggplot2', version = '3.5.1'); \
    remotes::install_version('ggrepel', version = '0.9.6'); \
    remotes::install_version('dplyr', version = '1.1.4'); \
    remotes::install_version('tidyr', version = '1.3.1'); \
    remotes::install_version('tibble', version = '3.2.1'); \
    remotes::install_version('purrr', version = '1.0.2'); \
    remotes::install_version('stringr', version = '1.5.1'); \
    remotes::install_version('forcats', version = '1.0.0'); \
    remotes::install_version('tidyverse', version = '2.0.0'); \
"

# Data handling packages
RUN Rscript -e " \
    remotes::install_version('data.table', version = '1.16.4'); \
    remotes::install_version('reshape2', version = '1.4.4'); \
    remotes::install_version('qs', version = '0.27.2'); \
"

# qs2 — faster qdata format for large object disk caching (successor to qs).
# Installed separately to preserve layer cache for the block above.
RUN Rscript -e "install.packages('qs2', repos = 'https://cloud.r-project.org')"

# Visualization packages
RUN Rscript -e " \
    remotes::install_version('viridis', version = '0.6.5'); \
    remotes::install_version('gridExtra', version = '2.3'); \
    remotes::install_version('scales', version = '1.3.0'); \
    remotes::install_version('cowplot', version = '1.1.3'); \
    remotes::install_version('egg', version = '0.4.5'); \
    remotes::install_version('ggpubr', version = '0.6.0'); \
    remotes::install_version('ggcorrplot', version = '0.1.4.1'); \
    remotes::install_version('ggplotify', version = '0.1.2'); \
    remotes::install_version('see', version = '0.9.0'); \
    remotes::install_version('svglite', version = '2.1.3'); \
    remotes::install_version('CMplot', version = '4.5.1'); \
    remotes::install_version('scattermore', version = '1.2'); \
    remotes::install_version('RColorBrewer', version = '1.1-3'); \
    remotes::install_version('VennDiagram', version = '1.8.2'); \
"

# Spatial and geographic packages
RUN Rscript -e " \
    remotes::install_version('sp', version = '2.1-4'); \
    remotes::install_version('terra', version = '1.8-5'); \
    remotes::install_version('raster', version = '3.6-30'); \
    remotes::install_version('geodata', version = '0.6-2'); \
    remotes::install_version('geosphere', version = '1.5-20'); \
    remotes::install_version('scatterpie', version = '0.2.4'); \
    remotes::install_version('ggnewscale', version = '0.5.0'); \
    remotes::install_version('ggspatial', version = '1.1.9'); \
"

# Shiny packages (legacy app)
RUN Rscript -e " \
    remotes::install_version('shinydashboard', version = '0.7.2'); \
    remotes::install_version('DT', version = '0.33'); \
    remotes::install_version('plotly', version = '4.10.4'); \
    remotes::install_version('htmlwidgets', version = '1.6.4'); \
"

# Shiny app (golem-based rewrite) dependencies
RUN Rscript -e " \
    remotes::install_version('bslib',     version = '0.9.0'); \
    remotes::install_version('bsicons',   version = '0.1.2'); \
    remotes::install_version('golem',     version = '0.5.1'); \
    remotes::install_version('yaml',      version = '2.3.10'); \
    remotes::install_version('jsonlite',  version = '1.8.9'); \
    remotes::install_version('base64enc', version = '0.1-3'); \
    remotes::install_version('cachem',    version = '1.1.0'); \
    remotes::install_version('config',    version = '0.3.2'); \
    remotes::install_version('rlang',     version = '1.1.4'); \
    remotes::install_version('thematic',  version = '0.1.5'); \
    remotes::install_version('processx',  version = '3.8.4'); \
    remotes::install_version('shinyjs',     version = '2.1.0'); \
    remotes::install_version('shinyFiles',  version = '0.9.3'); \
"

# Install golem Shiny app as R package
COPY scripts/clinego.app /tmp/clinego.app
RUN Rscript -e "remotes::install_local('/tmp/clinego.app', dependencies = FALSE)" \
  && rm -rf /tmp/clinego.app

# topr - CRITICAL: version >= 2.0.0 required for custom (non-human) genome builds
RUN Rscript -e "remotes::install_version('topr', version = '2.0.2')"

# BiocManager for Bioconductor packages
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install(version = '3.22', ask = FALSE)"

# Every BiocManager::install() below is followed by a requireNamespace() check that
# stops the layer. BiocManager::install() WARNS on a failed install and returns
# normally, so a bare RUN exits 0, BuildKit caches the broken layer as a success,
# and every later build reuses it — the image is then silently missing a package
# and only fails months later inside a pipeline rule. Seen twice on one image:
# a cached layer that never installed LEA, and another that never installed WGCNA.
# remotes::install_version() already throws, so only these lines need the guard.

# Bioconductor packages (pinned to 3.22 versions)
RUN Rscript -e "BiocManager::install('LEA', version = '3.22', ask = FALSE); if (!requireNamespace('LEA', quietly = TRUE)) stop('FAILED to install LEA')"
RUN Rscript -e "BiocManager::install('qvalue', version = '3.22', ask = FALSE); if (!requireNamespace('qvalue', quietly = TRUE)) stop('FAILED to install qvalue')"
RUN Rscript -e "BiocManager::install('GenomicRanges', version = '3.22', ask = FALSE); if (!requireNamespace('GenomicRanges', quietly = TRUE)) stop('FAILED to install GenomicRanges')"
RUN Rscript -e "BiocManager::install('WGCNA', version = '3.22', ask = FALSE); if (!requireNamespace('WGCNA', quietly = TRUE)) stop('FAILED to install WGCNA')"
RUN Rscript -e "BiocManager::install('clusterProfiler', version = '3.22', ask = FALSE); if (!requireNamespace('clusterProfiler', quietly = TRUE)) stop('FAILED to install clusterProfiler')"
RUN Rscript -e "BiocManager::install('AnnotationDbi', version = '3.22', ask = FALSE); if (!requireNamespace('AnnotationDbi', quietly = TRUE)) stop('FAILED to install AnnotationDbi')"
RUN Rscript -e "BiocManager::install('GO.db', version = '3.22', ask = FALSE); if (!requireNamespace('GO.db', quietly = TRUE)) stop('FAILED to install GO.db')"
RUN Rscript -e "BiocManager::install('enrichplot', version = '3.22', ask = FALSE); if (!requireNamespace('enrichplot', quietly = TRUE)) stop('FAILED to install enrichplot')"

# CRAN packages used with Bioconductor workflows
RUN Rscript -e " \
    remotes::install_version('vegan', version = '2.6-8'); \
    remotes::install_version('robust', version = '0.7-5'); \
    remotes::install_version('SpatialPack', version = '0.4-1'); \
    remotes::install_version('vcfR', version = '1.15.0'); \
    remotes::install_version('adegenet', version = '2.1.10'); \
    remotes::install_version('poppr', version = '2.9.6'); \
    remotes::install_version('ggraph', version = '2.2.1'); \
"

# Spatial eigenvector analysis (dbMEM/MEM) — preGEA varpart block + the SHARED
# spatial-vector artifact reused by RDA, varpart and Gradient Forest
# (docs/rda_research.md C4/A20). Deliberately a SEPARATE RUN layer: adespatial
# pulls a large dependency chain (ade4, adegraphics, adephylo, sp, spdep) and a
# failure here should not invalidate the vegan/robust layer above. Installed
# AFTER vegan 2.6-8 with upgrade='never' so it resolves against the PINNED
# vegan rather than silently pulling a newer one.
# System deps already present: libgdal-dev/libgeos-dev/libproj-dev/libudunits2-dev.
RUN Rscript -e "remotes::install_version('adespatial', version = '0.3-29', upgrade = 'never')"

# Haplotype analysis packages
RUN Rscript -e "remotes::install_version('crosshap', version = '1.4.0')"

# GAPIT3 association models (GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK)
RUN Rscript -e "BiocManager::install('multtest', version = '3.22', ask = FALSE); if (!requireNamespace('multtest', quietly = TRUE)) stop('FAILED to install multtest')"
RUN Rscript -e " \
    remotes::install_version('EMMREML', version = '3.1'); \
    remotes::install_version('bigmemory', version = '4.6.4'); \
    remotes::install_version('gplots', version = '3.1.3.1'); \
    remotes::install_version('scatterplot3d', version = '0.3-44'); \
    remotes::install_version('snowfall', version = '1.84-6.3'); \
"
RUN Rscript -e "remotes::install_github('jiabowang/GAPIT@GAPIT3.5')"

# PopLDdecay v3.43 (pinned) — fast LD decay computation from VCF
RUN git clone --branch v3.43 --depth 1 https://github.com/hewm2008/PopLDdecay.git /tmp/PopLDdecay \
    && cd /tmp/PopLDdecay/src && make \
    && cp /tmp/PopLDdecay/bin/PopLDdecay /usr/local/bin/PopLDdecay \
    && rm -rf /tmp/PopLDdecay

# gradientForest from r-forge (version 0.1-37)
# Requires patching for R 4.5 (Calloc/Free -> R_Calloc/R_Free)
RUN git clone --depth 1 https://github.com/r-forge/gradientforest.git /tmp/gf && \
    cd /tmp/gf/pkg/extendedForest/src && \
    sed -i 's/\bCalloc(/R_Calloc(/g; s/\bFree(/R_Free(/g' *.c && \
    cd /tmp/gf/pkg/gradientForest/src && \
    sed -i 's/\bCalloc(/R_Calloc(/g; s/\bFree(/R_Free(/g' *.c 2>/dev/null || true && \
    Rscript -e "install.packages('/tmp/gf/pkg/extendedForest', repos = NULL, type = 'source')" && \
    Rscript -e "install.packages('/tmp/gf/pkg/gradientForest', repos = NULL, type = 'source')" && \
    rm -rf /tmp/gf

# Re-pin vegan to 2.6-8 as the LAST package-installation step. adespatial's
# and crosshap's OWN dependency resolvers each silently pulled vegan up to a
# newer CRAN-current release (2.7-3) despite upgrade='never' on their own
# install calls — 'never' governs whether an ALREADY-SATISFIED dependency is
# upgraded for THAT package's install, not whether a LATER, unrelated
# package's install pass re-resolves "old packages" globally (verified:
# crosshap's dependency install logged "Old packages: poppr, vcfR, vegan" and
# upgraded vegan again after an earlier re-pin placed right after adespatial).
# Placing this after every other package install (nothing installs anything
# afterward) is what actually makes the pin stick. Caught by the
# packageVersion('vegan') canary in the verification RUN below — that check
# is exactly what caught both silent upgrades during development.
RUN Rscript -e "remotes::install_version('vegan', version = '2.6-8', upgrade = 'never')"

# Verify EVERY package the pipeline actually loads.
#
# The per-line requireNamespace() guards above catch a BiocManager::install()
# that warned instead of failing, but only on the lines that remember to carry
# one, and this final RUN used to check six hand-picked packages, none of them
# Bioconductor. Two things slip through that: a BiocManager line added later
# without a guard, and a package that arrives only as a transitive dependency
# of a pinned package and stops arriving when that package's DESCRIPTION
# changes (ggdist reaches this image solely as an Imports of crosshap 1.4.0 --
# present today, asserted by nothing before this check).
#
# The list below is the set of packages that scripts/**/*.R library() or ::
# (comments stripped, base+recommended and test-only deps excluded).
# Regenerate it the same way when adding a script; a package used by a script
# and absent here is exactly the hole this check exists to close.
#
# clinego.app is in the list even though it is not a dependency but the thing
# this Dockerfile builds. It belongs here because the failure it catches already
# happened: the COPY + install_local pair above produced a cached 0 B layer, so
# the image shipped for weeks with no app package while this very RUN reported
# every package present. `run_app()` was dead; only the bind-mounted dev.R path
# worked, which is why nobody noticed. Asserting the package the build installs
# is the cheapest way to keep a poisoned cache entry from passing as a build.
RUN Rscript -e " \
    required <- c( \
        'AnnotationDbi', 'DT', 'GAPIT', 'GO.db', 'LEA', 'SpatialPack', \
        'VennDiagram', 'WGCNA', 'adegenet', 'adespatial', 'base64enc', \
        'bsicons', 'bslib', 'cachem', 'clinego.app', 'clusterProfiler', 'config', \
        'cowplot', 'crosshap', 'data.table', 'digest', 'dplyr', \
        'enrichplot', 'forcats', 'geodata', 'geosphere', 'ggcorrplot', \
        'ggdist', 'ggnewscale', 'ggplot2', 'ggpubr', 'ggrepel', \
        'ggspatial', 'golem', 'gradientForest', 'htmltools', \
        'htmlwidgets', 'jsonlite', 'magrittr', 'plotly', 'poppr', \
        'processx', 'purrr', 'qs', 'qs2', 'qvalue', \
        'robust', 'sass', 'scales', 'scattermore', 'scatterpie', \
        'shiny', 'shinyFiles', 'shinyjs', 'stringr', 'svglite', \
        'terra', 'tibble', 'tidyr', 'topr', 'vcfR', 'vegan', 'viridis', \
        'yaml' \
    ); \
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]; \
    if (length(missing)) stop('MISSING R packages in image: ', paste(missing, collapse = ', ')); \
    cat('All', length(required), 'required packages present.\n'); \
    stopifnot(packageVersion('topr') >= '2.0.0'); \
    stopifnot(packageVersion('ggplot2') >= '3.5.0'); \
    stopifnot(packageVersion('ggrepel') >= '0.9.0'); \
    stopifnot(requireNamespace('gradientForest', quietly = TRUE)); \
    stopifnot(requireNamespace('crosshap', quietly = TRUE)); \
    stopifnot(requireNamespace('GAPIT', quietly = TRUE)); \
    library(robust); \
    x <- matrix(rnorm(300), ncol = 3); \
    stopifnot(length(covRob(x, distance = TRUE, na.action = na.omit, estim = 'pairwiseGK')\$dist) == 100); \
    stopifnot(packageVersion('vegan') == '2.6.8'); \
    library(adespatial); \
    set.seed(1); xy <- cbind(runif(12), runif(12)); \
    .d <- dist(xy); .th <- max(vegan::spantree(.d)\$dist); \
    .m <- dbmem(.d, thresh = .th, MEM.autocor = 'positive', silent = TRUE); \
    stopifnot(ncol(as.data.frame(.m)) >= 1); \
    cat('All package version checks passed.\n'); \
"
