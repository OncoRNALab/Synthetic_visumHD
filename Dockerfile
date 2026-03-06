# Visium-HD-like spatial FASTQ simulation: sCCIgen + ART + Python
# Build: docker build -t stdata_sim .
FROM rocker/r-ver:4.2.0

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    gzip \
    ca-certificates \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R dependencies and sCCIgen (optional; pipeline has R fallback)
RUN R -e "install.packages(c('devtools', 'Matrix', 'rhdf5'), repos='https://cloud.r-project.org')" || true
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install('rhdf5', ask=FALSE)" || true
RUN R -e "devtools::install_github('songxiaoyu/sCCIgen', dependencies=TRUE)" || true

# Install ART if available (optional; pipeline uses Python fallback)
# USER can mount ART binary or install in derived image
# RUN apt-get install -y art-nextgen-simulation-tools || true

WORKDIR /data
COPY bin/ /data/bin/
COPY scripts/ /data/scripts/
ENV PATH="/data/bin:${PATH}"

# Nextflow runs processes in containers; projectDir is mounted
ENTRYPOINT ["/bin/bash"]
CMD ["-c", "echo 'Use with Nextflow -profile docker'"]
