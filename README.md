# 🧬 Methyline

> A modular and reproducible Nextflow pipeline for Enzymatic Methyl-seq (EM-seq) data analysis

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Overview

**Methyline** is a modular, reproducible bioinformatics pipeline built with [Nextflow DSL2](https://www.nextflow.io/) for the end-to-end analysis of Enzymatic Methyl-seq (EM-seq) data. It is designed to be applicable to both clinical and research settings, prioritising reproducibility, portability, and ease of use.

The pipeline covers the full analysis workflow, from raw read quality control through to functional annotation of differentially methylated regions (DMRs) and result visualisation.

> ⚠️ **This pipeline is under active development.** Tool selection and module implementation are ongoing. The current repository reflects the architectural design; functional modules will be progressively integrated.

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────┐
│                        METHYLINE                        │
│          EM-seq Analysis Pipeline (Nextflow DSL2)       │
└─────────────────────────────────────────────────────────┘

  Raw reads (FASTQ)
       │
       ▼
┌─────────────┐
│  1. QC &    │  Quality control and adapter/quality trimming
│  Trimming   │
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  2. Align-  │  Bisulfite-aware alignment to reference genome
│   ment      │
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  3. Methyl- │  CpG methylation extraction and quantification
│  Extraction │
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  4. DMR     │  Identification of differentially methylated
│  Analysis   │  regions across sample groups
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  5. Funct.  │  Genomic annotation and functional enrichment
│  Annotation │  of DMRs
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  6. Visual- │  Automated report generation and result
│  isation    │  visualisation
└─────────────┘
```

---

## Features

- **Modular design**: each analysis step is implemented as an independent Nextflow module, allowing flexible execution and easy substitution of tools
- **Reproducible**: full containerisation support via Docker and Singularity
- **Portable**: runs locally, on HPC clusters (SLURM, LSF), or cloud environments (AWS, Google Cloud)
- **EM-seq optimised**: tool selection is specifically evaluated for compatibility with EM-seq data characteristics
- **Automated reporting**: integrated generation of quality and methylation summary reports

---

## Pipeline Modules

| Step | Module | Description | Status |
|------|--------|-------------|--------|
| 1 | QC & Trimming | Raw read quality assessment and adapter trimming | 🔧 In development |
| 2 | Alignment | Bisulfite-aware alignment to reference genome | 🔧 In development |
| 3 | Methylation Extraction | CpG methylation calling and quantification | 🔧 In development |
| 4 | DMR Analysis | Differential methylation analysis across conditions | 🔧 In development |
| 5 | Functional Annotation | Annotation of DMRs with genomic features and pathway enrichment | 🔧 In development |
| 6 | Visualisation | QC plots, methylation profiles and summary HTML report | 🔧 In development |

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/) (recommended for reproducibility)
- Java >= 11

---

## Quick Start

> ⚠️ Execution instructions will be updated as modules are implemented.

```bash
# Clone the repository
git clone https://github.com/<your-username>/methyline.git
cd methyline

# Run the pipeline (Docker profile)
nextflow run main.nf \
  --input samplesheet.csv \
  --genome hg38 \
  --outdir results/ \
  -profile docker
```

---

## Input

The pipeline accepts a **samplesheet** in CSV format specifying sample identifiers, paths to FASTQ files, and experimental group metadata. Full samplesheet specification will be documented upon module completion.

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier |
| `fastq_1` | Path to R1 FASTQ file |
| `fastq_2` | Path to R2 FASTQ file (paired-end) |
| `condition` | Experimental group (e.g. `case`, `control`) |

---

## Output

The pipeline generates a structured output directory organised by module:

```
results/
├── qc/                  # FastQC / MultiQC reports
├── trimming/            # Trimmed reads
├── alignment/           # BAM files and alignment statistics
├── methylation/         # CpG methylation tables (bedGraph / bismark coverage)
├── dmr/                 # DMR tables and summary statistics
├── annotation/          # Annotated DMR files
└── report/              # Final HTML summary report
```

---

## Configuration

Pipeline behaviour can be customised via `nextflow.config` and parameter files. Key parameters (reference genome, tool-specific settings, resource allocation) will be documented as development progresses.

---

## Citation

If you use Methyline in your research, please cite:

> [Citation will be added upon publication]

---

## Author

Developed as part of a Master's thesis in Bioinformatics.  
Contact: negido@uoc.edu
