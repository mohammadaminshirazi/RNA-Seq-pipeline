# ==============================================================================
# Dockerfile for RNA-seq Pipeline
# ==============================================================================
# Build: docker build -t rnaseq-pipeline .
# Run: docker run -v /path/to/data:/data rnaseq-pipeline
# ==============================================================================

FROM condaforge/mambaforge:latest

LABEL maintainer="Mohammad Amin Shirazi"
LABEL description="RNA-seq Differential Expression Pipeline"
LABEL version="1.0.0"

# Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean -afy

# Activate environment by default
SHELL ["conda", "run", "-n", "rnaseq-pipeline", "/bin/bash", "-c"]

# Copy pipeline files
WORKDIR /pipeline
COPY . /pipeline/

# Make scripts executable
RUN chmod +x scripts/*.sh scripts/*.py

# Set entrypoint
ENTRYPOINT ["conda", "run", "-n", "rnaseq-pipeline"]
CMD ["snakemake", "--help"]
