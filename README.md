# 🧬 Reproducible RNA-seq Differential Expression Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![R 4.0+](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

A production-ready, reproducible bioinformatics pipeline for bulk RNA-seq differential gene expression analysis. Built with Snakemake for workflow management, featuring comprehensive quality control, multiple alignment options, and publication-ready visualizations.

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Pipeline Architecture](#pipeline-architecture)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output Description](#output-description)
- [Example Analysis](#example-analysis)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

## 🎯 Overview

This pipeline performs end-to-end analysis of bulk RNA-seq data, from raw FASTQ files to differential expression results and publication-ready visualizations. It is designed to be:

- **Reproducible**: Conda environments and containerization ensure consistent results
- **Scalable**: Efficiently handles datasets from a few samples to hundreds
- **Flexible**: Configurable parameters for different experimental designs
- **Well-documented**: Comprehensive documentation and example datasets

### Biological Applications

- Cancer vs. Normal tissue comparisons
- Treatment vs. Control experiments
- Developmental stage comparisons
- Disease progression studies
- Drug response analysis

## ✨ Features

### Quality Control
- **FastQC**: Per-sample quality metrics
- **MultiQC**: Aggregated quality reports
- **Adapter trimming**: Automatic adapter detection and removal (Trimmomatic/fastp)

### Alignment & Quantification
- **STAR**: Splice-aware alignment (recommended for most applications)
- **HISAT2**: Memory-efficient alternative
- **featureCounts**: Gene-level quantification
- **Salmon**: Transcript-level quantification (optional)

### Differential Expression
- **DESeq2**: Primary DE analysis engine
- **edgeR**: Alternative DE method for comparison
- **Multiple testing correction**: Benjamini-Hochberg FDR

### Visualizations
- Volcano plots with customizable thresholds
- MA plots (log ratio vs. mean expression)
- Hierarchical clustering heatmaps
- PCA plots for sample relationships
- Gene expression boxplots
- Interactive HTML reports

### Functional Analysis
- Gene Ontology (GO) enrichment
- KEGG pathway analysis
- Gene Set Enrichment Analysis (GSEA)

## 🏗️ Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        RNA-seq Analysis Pipeline                        │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Raw FASTQ │────▶│     QC      │────▶│   Trimming  │────▶│  Alignment  │
│    Files    │     │  (FastQC)   │     │(Trimmomatic)│     │   (STAR)    │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
                                                                   │
                                                                   ▼
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Results   │◀────│  DESeq2/    │◀────│   Count     │◀────│     BAM     │
│   & Plots   │     │   edgeR     │     │   Matrix    │     │    Files    │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  📊 Outputs: Volcano, MA, Heatmap, PCA, GO/KEGG, HTML Report           │
└─────────────────────────────────────────────────────────────────────────┘
```

## 🚀 Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/rnaseq-pipeline.git
cd rnaseq-pipeline

# 2. Create and activate conda environment
conda env create -f environment.yml
conda activate rnaseq-pipeline

# 3. Configure your analysis
cp config/config_template.yaml config/config.yaml
# Edit config.yaml with your sample information

# 4. Run the pipeline (dry-run first)
snakemake --dry-run --cores 4

# 5. Execute the pipeline
snakemake --cores 8 --use-conda
```

## 📦 Installation

### Prerequisites

- Linux/macOS operating system
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- At least 16GB RAM (32GB recommended for large datasets)
- ~50GB disk space for reference genome and indices

### Step-by-step Installation

```bash
# Option 1: Using Conda
conda env create -f environment.yml
conda activate rnaseq-pipeline

# Option 2: Using Mamba (faster)
mamba env create -f environment.yml
mamba activate rnaseq-pipeline

# Option 3: Docker (most reproducible)
docker pull YOUR_USERNAME/rnaseq-pipeline:latest
docker run -v /path/to/data:/data rnaseq-pipeline
```

### Reference Genome Setup

```bash
# Download reference genome and annotation (example: human GRCh38)
bash scripts/download_references.sh --genome hg38 --output resources/

# Or for mouse (GRCm39)
bash scripts/download_references.sh --genome mm39 --output resources/
```

## ⚙️ Configuration

### Sample Sheet Format

Create a `samples.tsv` file with the following format:

```
sample_id	condition	replicate	fastq_1	fastq_2
tumor_1	tumor	1	data/tumor_1_R1.fastq.gz	data/tumor_1_R2.fastq.gz
tumor_2	tumor	2	data/tumor_2_R1.fastq.gz	data/tumor_2_R2.fastq.gz
tumor_3	tumor	3	data/tumor_3_R1.fastq.gz	data/tumor_3_R2.fastq.gz
normal_1	normal	1	data/normal_1_R1.fastq.gz	data/normal_1_R2.fastq.gz
normal_2	normal	2	data/normal_2_R1.fastq.gz	data/normal_2_R2.fastq.gz
normal_3	normal	3	data/normal_3_R1.fastq.gz	data/normal_3_R2.fastq.gz
```

### Configuration File (config.yaml)

```yaml
# See config/config.yaml for full documentation
samples: "config/samples.tsv"
reference:
  genome: "resources/genome.fa"
  annotation: "resources/annotation.gtf"
  
analysis:
  contrast: ["tumor", "normal"]
  fdr_threshold: 0.05
  log2fc_threshold: 1.0
```

## 📖 Usage

### Basic Usage

```bash
# Full pipeline execution
snakemake --cores 8 --use-conda

# Generate specific outputs
snakemake results/deseq2/differential_expression.csv --cores 4

# Generate report
snakemake --report report.html
```

### Advanced Options

```bash
# Dry-run to see what will be executed
snakemake --dry-run

# Run with cluster submission (SLURM)
snakemake --cluster "sbatch -p normal -t {resources.time}" --jobs 100

# Run with specific config
snakemake --configfile config/custom_config.yaml --cores 8

# Force re-run of specific rules
snakemake --forcerun deseq2_analysis --cores 8
```

## 📁 Output Description

```
results/
├── qc/
│   ├── fastqc/                    # Individual FastQC reports
│   ├── multiqc/                   # Aggregated QC report
│   └── multiqc_report.html        # Interactive QC dashboard
├── alignment/
│   ├── {sample}.sorted.bam        # Aligned reads
│   ├── {sample}.sorted.bam.bai    # BAM index
│   └── alignment_stats.txt        # Alignment statistics
├── counts/
│   ├── raw_counts.csv             # Raw gene counts
│   ├── normalized_counts.csv      # DESeq2 normalized counts
│   └── tpm_counts.csv             # TPM normalized counts
├── deseq2/
│   ├── differential_expression.csv # Full DE results
│   ├── significant_genes.csv      # Filtered significant genes
│   ├── deseq2_results.rds         # R object for further analysis
│   └── analysis_summary.txt       # Summary statistics
└── visualizations/
    ├── volcano_plot.pdf           # Volcano plot
    ├── volcano_plot.png           # PNG version
    ├── ma_plot.pdf                # MA plot
    ├── heatmap_top_genes.pdf      # Heatmap of top DE genes
    ├── pca_plot.pdf               # PCA plot
    ├── sample_correlation.pdf     # Sample correlation heatmap
    └── interactive_report.html    # Interactive HTML report
```

## 📊 Example Analysis

### Using Example Data

```bash
# Download example dataset (breast cancer vs normal)
bash scripts/download_example_data.sh

# Run pipeline with example config
snakemake --configfile config/example_config.yaml --cores 4
```

### Expected Results

After successful completion, you should see:
- ~15,000-20,000 expressed genes
- ~500-3000 differentially expressed genes (FDR < 0.05, |log2FC| > 1)
- Publication-ready figures in `results/visualizations/`

## 🔧 Troubleshooting

### Common Issues

1. **Memory errors during STAR alignment**
   ```bash
   # Increase memory or use HISAT2
   snakemake --config aligner=hisat2 --cores 8
   ```

2. **Missing dependencies**
   ```bash
   # Recreate environment
   conda env remove -n rnaseq-pipeline
   conda env create -f environment.yml
   ```

3. **File not found errors**
   - Check sample sheet paths are correct
   - Ensure FASTQ files exist and are readable

### Getting Help

- Open an issue on GitHub
- Check existing issues for solutions
- Contact: your.email@example.com

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📚 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{rnaseq_pipeline,
  author = {Your Name},
  title = {Reproducible RNA-seq Differential Expression Pipeline},
  year = {2025},
  url = {https://github.com/YOUR_USERNAME/rnaseq-pipeline}
}
```

Also cite the tools used:
- DESeq2: Love MI, Huber W, Anders S (2014). Genome Biology 15:550
- STAR: Dobin A et al. (2013). Bioinformatics 29:15-21
- Snakemake: Mölder F et al. (2021). F1000Research 10:33

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Made with ❤️ for the bioinformatics community**

*Last updated: February 2025*
