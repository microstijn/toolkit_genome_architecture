# Genome Architecture & Taxonomy Toolkit

This toolkit provides a suite of modular, command-line driven scripts for analyzing prokaryotic genomes. The primary workflow extracts taxonomic and assembly metadata from NCBI datasets, retrieves associated environmental data from the OmniMicrobe API, and calculates detailed gene overlap and genome architecture metrics from GFF files.

The entire process is managed by a central pipeline controller (`run_pipeline.jl`) for automated, end-to-end analysis.

## Features

* **Automated Pipeline**: A master script (`run_pipeline.jl`) runs the entire workflow, from raw data to final analysis, in the correct order.
* **Comprehensive Analysis**:
    * Extracts full taxonomic lineage (superkingdom to species).
    * Maps organisms to known environments using the OmniMicrobe database.
    * Calculates unidirectional, convergent, and divergent gene overlaps.

## Workflow

The master script `run_pipeline.jl` manages this flow automatically.

---
## Data Sources

All data processed by this toolkit is sourced from the National Center for Biotechnology Information (NCBI).

### Genome Assemblies and Metadata
The required genome assembly and metadata files can be downloaded using the NCBI `datasets` command-line tool. This process retrieves GFF3 files for genome annotation and the accompanying JSONL files (`assembly_data_report.jsonl`) that are processed by the scripts.

For example, to download all available reference genomes for **Bacteria (taxon: 2)** and **Archaea (taxon: 2157)**, you can use commands like the following:
```bash
# Download dehydrated packages for Bacteria and Archaea
datasets download genome taxon 2 --assembly-source refseq --reference --include gff3,gtf,seq-report --dehydrated --filename bacteria_reference.zip
datasets download genome taxon 2157 --assembly-source refseq --reference --include gff3,gtf,seq-report --dehydrated --filename archaea_reference.zip

# Rehydrate the packages to retrieve the data files
datasets rehydrate --directory bacteria_reference/
datasets rehydrate --directory archaea_reference/
```

### Taxonomy Database
The taxonomic lineage information is derived from the NCBI Taxonomy database dump files (`nodes.dmp` and `names.dmp`). These files can be obtained from the following FTP address:

`ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`

---

## Installation

### Prerequisites

* **Julia**: Version 1.6 or higher.

### Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/microstijn/toolkit_genome_architecture.git
    cd toolkit_genome_architecture/src
    ```

2.  **Install Julia Dependencies:**
    This project uses a local Julia environment defined by `Project.toml` and `Manifest.toml`. To install all required packages with their correct versions, start Julia in the `src` directory and run the following commands:
    ```julia
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    ```

---

## Usage

### Recommended Method: The Master Pipeline

The easiest and recommended way to use the toolkit is via the `run_pipeline.jl` script. It handles all steps, file management, and error checking.

Run it from your terminal within the `src` directory, providing paths to your data.

**Example Command:**

```bash
julia run_pipeline.jl \
    --jsonl-files /path/to/your/data/bacteria_report.jsonl /path/to/your/data/archaea_report.jsonl \
    --gff-dir /path/to/your/gff_files/ \
    --taxdump-dir /path/to/your/taxdump/ \
    --output-dir ./pipeline_output
```

**Arguments:**

* `--jsonl-files`: Path(s) to your NCBI `assembly_data_report.jsonl` files.
* `--gff-dir`: Path to the directory containing all your `.gff` genome files.
* `--taxdump-dir`: Path to the directory containing the NCBI `nodes.dmp` and `names.dmp` files.
* `--output-dir`: The directory where all results will be saved. It will be created if it doesn't exist.

### Running Individual Scripts

You can also run each script individually. Use the `--help` flag to see all available options for each script (e.g., `julia taxid_from_ncbi_JSON.jl --help`).

---

## Scripts Overview

#### `run_pipeline.jl`

* **Description**: The master controller script that orchestrates the entire workflow. This is the recommended entry point.
* **Inputs**: Command-line paths to all raw data directories.
* **Outputs**: A populated output directory containing all generated analysis files.

#### `taxid_from_ncbi_JSON.jl`

* **Description**: Parses NCBI `assembly_data_report.jsonl` files to extract assembly metadata (accession, completeness, GC content) and organism info. It then uses an NCBI taxdump to map each entry to its full taxonomic lineage.
* **Inputs**: NCBI JSONL file(s) and the `nodes.dmp`/`names.dmp` taxdump files.
* **Outputs**: A CSV file (e.g., `bacteria_report_TaxId.csv`) containing the combined metadata and taxonomy.

#### `environment_from_omnicrobe_by_taxid.jl`

* **Description**: Takes the taxonomy CSV from the previous step and queries the OmniMicrobe API to find known environmental niches for each taxon ID. It uses concurrent (asynchronous) requests to speed up the process.
* **Inputs**: The `...TaxId.csv` file.
* **Outputs**: A new CSV file (e.g., `...TaxId_Omni.csv`) with added columns for `environments` and `obtId`.

#### `calc_gene_overlap_on_GFF.jl`

* **Description**: Analyzes a directory of GFF3 files to calculate detailed genome architecture metrics. It processes each contig to measure gene density, spacing, and overlaps (unidirectional, convergent, and divergent).
* **Inputs**: A directory containing GFF3 files.
* **Outputs**: A single CSV file (`genome_architecture_metrics.csv`) summarizing the structural features of every contig in every input genome.

---

## Author & Revision

* **Original Author**: SHP (2022-2025)
* **Last Revised**: August 7, 2025
* **License**: This project is licensed under the MIT License.
