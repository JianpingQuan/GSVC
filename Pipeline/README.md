# GSVC bioinformatics pipeline

This directory is reserved for the **upstream metagenomic pipeline** used to build the Global Swine Virome Catalogue (GSVC):

1. Read QC and host depletion (`fastp`, `Bowtie2`)
2. Assembly (`MEGAHIT`)
3. MAG recovery (`MetaBAT2`, `CheckM`, `dRep`, `GTDB-Tk`)
4. Viral identification (`VirSorter2`, `VIBRANT`, `viralVerify`, `geNomad`, `CheckV`)
5. vOTU clustering (`MMseqs2`, 95% ANI / 85% AF)
6. Abundance profiling (`CoverM` → Count / TPM matrices)
7. Host prediction (`iPHoP`), CRISPR spacer linking, functional annotation

## Contents

| Path | Description |
|------|-------------|
| `modules/` | Shell wrappers for QC, host filtering, assembly, binning, viral identification, CoverM |
| `sbatch/` | SLURM batch templates (adapt paths before HPC submission) |
| `script/` | Python/R utilities (preprocessing, CRISPR merging, lifestyle integration, CoverM merge) |
| `pipeline.txt` | End-to-end command log from catalogue construction |
| `bioproject.list.txt` | BioProject accessions processed |

> **Note:** Paths in `pipeline.txt` and sbatch files point to the authors' HPC environment. Replace `/work/home/...` with your local paths, or wrap modules in Snakemake/Nextflow for portability.

For downstream figure reproduction, use:

- **`Scripts/`** — statistical analyses and plotting
- **`Pre-processed_Files/`** — harmonized metadata and summary tables
- **`large_data/`** — Zenodo-hosted Count/TPM matrices

## Reproducing the catalogue

See the root `README.md` and the manuscript *Materials and Methods*.
