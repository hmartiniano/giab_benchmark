# GIAB Data Download Workflow

This Snakemake workflow automates the acquisition of Genome in a Bottle (GIAB) reference materials, including FASTA genomes, VCF truth sets, and associated index files (.fai, .gzi, .tbi, .bed). Designed for reproducibility and scalability in bioinformatics pipelines.

## Purpose

Facilitates the standardized download and organization of GIAB benchmark data, crucial for variant calling pipeline development and validation. It abstracts complex FTP/HTTPS data retrieval into a declarative workflow.

## Key Features

*   **Dynamic Target Generation**: `ALL_TARGETS` are dynamically constructed from `config.yaml`, ensuring all specified data components (reference, VCF, TBI, BED, FASTQ) are targeted.
*   **Modular Download Rules**: Dedicated Snakemake rules (`download_reference`, `download_giab_vcf`, etc.) handle specific data types, promoting reusability and clarity.
*   **Config-Driven**: All data sources, sample definitions (including `latest` release URLs), and output paths are managed via `config.yaml`, enabling easy customization without modifying the `Snakefile`.
*   **Index File Integration**: Automatically retrieves `.fai` and `.gzi` for reference genomes, and `.tbi` for VCFs.

## Data Sources

Data is sourced primarily from NCBI's GIAB FTP/HTTPS repositories and NIST's opendata platform. URLs are explicitly defined per sample in `config.yaml`.

## Dependencies

*   [Snakemake](https://snakemake.readthedocs.io/) (workflow management)
*   `wget` (data retrieval)

## Usage

Execute the workflow from the project root:

```bash
snakemake --cores <N>
```
Replace `<N>` with the number of available CPU cores for parallel downloads.

## Configuration (`config.yaml`)

The `config.yaml` is the central configuration file. Define `references` (including `url` and `filename`), and `samples`. Each `sample` entry requires:
*   `id`: Sample identifier (e.g., NA24385).
*   `giab_version`: GIAB benchmark version (e.g., v4.2.1).
*   `reference_build`: Reference genome build (e.g., GRCh38).
*   `reference_key`: Maps to an entry in the `references` section.
*   `vcf_url`, `tbi_url`, `bed_url`: Direct URLs to the respective GIAB truth set files. These should point to `latest` release directories for dynamic updates.
*   `fastq_urls`: Direct URLs to FASTQ files.

Commented-out sections in `config.yaml` demonstrate how to include additional samples (e.g., HG001, Chinese Trio) by uncommenting and providing valid URLs.

## Output Structure

Downloaded data is organized under the `data/` directory (configurable via `data_dir`), with subdirectories for `reference/`, `giab/`, and `fastq/` structured by sample name, GIAB version, and reference build.
