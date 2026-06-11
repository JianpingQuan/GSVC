# Construction of the Global Swine Virome Catalogue (GSVC) and metagenome-assembled genomes of the pig gut microbiome

This directory contains scripts related to the manuscript **"A metagenomic catalogue of the global swine gut virome uncovers ecological principles governing phage-bacteria interactions"** (in preparation).

The repository provides the computational workflow used to build the **Global Swine Virome Catalogue (GSVC)** — **224,605 vOTUs** from **5,064** pig gut metagenomes across **25 countries** — together with a companion set of **13,862** dereplicated bacterial metagenome-assembled genomes (MAGs) for CRISPR-based host–virus linking, downstream ecological analyses, and figure reproduction.

**Before running, you must ensure that all required software and databases are installed successfully.**

---

## INSTALLATION

Create two directories **`bin`** and **`Database`** in your home directory (or another location of your choice), and add `~/bin` to your `PATH`.

### Software installation

Install each tool following its official manual. Versions listed below are those used in the paper; newer versions are generally acceptable unless noted otherwise. Some tools bundle dependencies (e.g., Samtools with metaWRAP), so you do not need to install every dependency separately.

| Software | Version (paper) | Availability |
|----------|-----------------|--------------|
| fastp | v1.0 | https://github.com/OpenGene/fastp |
| Bowtie 2 | v2.5.2 | https://github.com/BenLangmead/bowtie2 |
| Samtools | v1.10+ | https://github.com/samtools/samtools |
| MEGAHIT | v1.2.9 | https://github.com/voutcn/megahit |
| MetaBAT2 | v2.15 | https://bitbucket.org/berkeleylab/metabat |
| CheckM | v1.2.4 | https://ecogenomics.github.io/CheckM/ |
| metaWRAP | v1.1.1 | https://github.com/bxlab/metaWRAP |
| dRep | v3.6.2 | https://github.com/MrOlm/drep |
| GTDB-Tk | release 226 | https://ecogenomics.github.io/GTDBTk/ |
| VirSorter2 | v2.2.3 | https://github.com/jiarong/VirSorter2 |
| VIBRANT | v1.2.1 | https://github.com/AnantharamanLab/VIBRANT |
| viralVerify | v1.1 | https://github.com/abai5214/viralVerify |
| geNomad | v1.5.0 | https://portal.nersc.gov/genomad/ |
| CheckV | v1.0.1 | https://bitbucket.org/berkeleylab/checkv |
| MMseqs2 | v18.8cc5c | https://github.com/soedinglab/MMseqs2 |
| CoverM | v0.6.1 | https://github.com/wwood/CoverM |
| Prodigal | v2.6.3 | https://github.com/hyattpd/Prodigal |
| iPHoP | v1.3.3 | https://bitbucket.org/berkeleylab/iphop |
| CRT | v1.2-CLI | https://github.com/davidw/devis/CRT |
| Piler-CR | v1.06 | https://www.drive5.com/pilercr/ |
| BLAST+ | v2.12.0+ | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| BACPHLIP | v0.9.6 | https://github.com/adamhockenberry/bacphlip |
| PhaTYP | v2.0 | https://github.com/gumeng999/PhaTYP |
| PhaBOX2 (PhaGCN) | v2.1.13 | https://github.com/KennthShang/PhaBOX |
| eggNOG-mapper | v2.1.12 | http://eggnog-mapper.embl.de/ |
| DIAMOND | v0.9.14.115 | https://github.com/bbuchfink/diamond |
| AMRFinderPlus | v3.11.2 | https://github.com/ncbi/amr |
| vContact3 | v3.1.6 | https://bitbucket.org/MAVERICLab/vcontact3 |
| seqkit | v2.13.0 | https://github.com/shenwei356/seqkit |
| R | ≥ 4.2 | https://www.r-project.org/ |
| Python | ≥ 3.9 | https://www.python.org/ |

**Note:** Make all required commands available in `~/bin` or in your system `PATH`. SLURM job templates in `Pipeline/sbatch/` assume a Linux HPC environment; adapt paths in `Pipeline/pipeline.txt` and `Pipeline/modules/` before submission. Bowtie2, Samtools, and MEGAHIT are also available through metaWRAP.

---

### Database installation

Store reference databases under **`~/Database`** (or set equivalent environment variables such as `GTDBTK_DATA_PATH` and geNomad DB paths).

| Database | Version / release | Description | Availability |
|----------|-------------------|-------------|--------------|
| Pig reference genome (*Sus scrofa*) | Sscrofa11.1 (or equivalent) | Host read removal | http://asia.ensembl.org/Sus_scrofa/Info/Index |
| CheckM / CheckV reference | bundled with tool | MAG / viral quality assessment | Installed with CheckM / CheckV |
| GTDB | release 226 | GTDB-Tk taxonomic assignment | https://gtdb.ecogenomic.org/downloads |
| geNomad database | v1.5 | Viral identification & taxonomy | Installed with geNomad (`genomad download-db`) |
| iPHoP database | Jun_2025_pub_rw | Host prediction | https://bitbucket.org/berkeleylab/iphop |
| eggNOG | EggNOG 5.0 | Functional orthology annotation | http://eggnog5.embl.de/#/app/downloads |
| SARG | v3.1 | Antibiotic resistance genes | https://smile.hku.hk/SARGs |
| VFDB | January 2026 release | Virulence factors | http://www.mgc.ac.cn/VFs/ |
| ICTV Virus Metadata Resource | MSL37.3 (v20-190822) | Viral novelty benchmarking | https://ictv.global/vmr |
| IMG/VR | v4 | Viral reference genomes | https://img.jgi.doe.gov/vr/ |
| VIBRANT databases | v1.2.1 | AMG and viral gene annotation | Installed with VIBRANT |
| VirSorter2 databases | v2.2.3 | Viral identification | Installed with VirSorter2 |

**Note:** Database versions reflect those used in the paper. Most public databases are updated regularly; results may differ slightly when using newer releases.

---

## OVERVIEW OF PIPELINE

Metagenomic analysis scripts are placed in the **`Pipeline`** directory. The workflow has two major modules — **(I) construction of the GSVC viral catalogue** and **(II) reconstruction of bacterial MAGs for CRISPR host linking** — that share the same preprocessing and assembly steps.

Command-level examples and batch submission logs are recorded in **`Pipeline/pipeline.txt`**. Shell wrappers are in **`Pipeline/modules/`**, SLURM templates in **`Pipeline/sbatch/`**, and helper scripts in **`Pipeline/script/`**.

---

### Shared steps (Modules I–II)

Steps before viral mining and MAG binning are identical for all samples.

| Part | Script / module | Description |
|------|-----------------|-------------|
| **Part 1** | `modules/00_raw_qc.sh` · `sbatch/00_raw_qc.sbatch` | **Metagenomic pre-processing:** adapter trimming and quality filtering with **fastp**; technical replicates merged when needed. |
| **Part 2** | `modules/01_host_filter.sh` · `modules/01_host_filter_pool.sh` | **Host (pig) read removal:** unmapped reads retained with **Bowtie2** against the *Sus scrofa* reference genome. |
| **Part 3** | MEGAHIT (via pipeline) | **Metagenomic assembly:** per-sample **de novo** assembly (minimum contig length 1,000 bp). |

---

### Construction of metagenome-assembled genomes (Module II)

MAGs support CRISPR spacer extraction and direct virus–host inference. High-quality MAGs were merged with published swine gut MAGs and dereplicated globally.

| Part | Script / module | Description |
|------|-----------------|-------------|
| **Part 4** | `modules/II-1_binning.sh` · `modules/metaWrap_binning.sh` | **Binning:** contigs binned with **MetaBAT2** (≥1,500 bp contigs). |
| **Part 5** | `modules/II-2_checkm.sh` | **Quality control:** bin completeness and contamination assessed with **CheckM** (lineage workflow). |
| **Part 6** | `modules/metaWrap_refinement.sh` · `script/II-3_filter_bins_median.py` · `script/II-4_rename_bins.sh` | **Refinement & filtering:** retain bins with completeness ≥90% and contamination <5%. |
| **Part 7** | `sbatch/II-5a_drep_dereplicate.sbatch` · `sbatch/II-5b_drep_dereplicate_combine.sbatch` | **Dereplication:** per-study and global dereplication with **dRep** (-comp 90 -con 5 -sa 0.99). |
| **Part 8** | `sbatch/II-6_gtdbtk.sbatch` | **Taxonomic classification:** **GTDB-Tk** classify workflow on dereplicated MAGs. |
| **Part 9** | `modules/coverM.sh` · `script/merge_coverm_multicol.py` | **MAG abundance:** read mapping with **CoverM**; TPM matrix generation. |
| **Part 10** | `script/identify_crispr.py` · `script/merge_crispr.py` | **CRISPR arrays:** spacer calling with **CRT** and **Piler-CR** on MAGs; spacers aligned to vOTU representatives with **BLASTN** (≥95% identity, ≤2 mismatches). |

**Output:** 13,862 non-redundant bacterial genomes used as the CRISPR spacer source database.

---

### Construction of the Global Swine Virome Catalogue (Module I)

Viral contigs were recovered from per-sample assemblies, quality-filtered, clustered into vOTUs, and functionally annotated.

| Part | Script / module | Description |
|------|-----------------|-------------|
| **Part 11** | `sbatch/virsorter2_viralPrediction.sbatch` | **VirSorter2** screening (≥1,500 bp contigs). |
| **Part 12** | `sbatch/vibrant_viralverify_viralPrediction.sbatch` | **VIBRANT** and **viralVerify** parallel screening. |
| **Part 13** | `sbatch/viral_identification.sbatch` | **Integrated viral identification:** composite scoring across three tools; contigs detected by all three tools require mean score >1.0, by two tools require mean score >4/3. |
| **Part 14** | geNomad · CheckV | **Quality filtering:** geNomad viral score >0.9, length ≥5 kb; CheckV completeness stratification (HQ ≥90%, MQ 50–90%) and provirus boundary refinement. |
| **Part 15** | MMseqs2 easy-cluster | **vOTU clustering:** ≥95% nucleotide identity, ≥85% alignment fraction (coverage mode 1); HQ representatives preferentially selected. |
| **Part 16** | geNomad (representatives) | **Final catalogue taxonomy:** representative sequences re-evaluated with geNomad (viral score ≥0.9). |
| **Part 17** | BACPHLIP · PhaTYP · `script/integrate_lifestyle_thresh095.py` | **Lifestyle inference:** temperate vs. virulent calls require concordance between **BACPHLIP** and **PhaTYP** (probability ≥0.95). |
| **Part 18** | iPHoP · `sbatch/iphop.sbatch` · `sbatch/iphop_hq.sbatch` | **Host prediction:** taxonomy assigned with **iPHoP** (Jun_2025_pub_rw reference). |
| **Part 19** | `modules/coverM.sh` · `script/universal_merge_coverm.py` | **vOTU abundance profiling:** **CoverM** mapping to GSVC representatives; global Count and TPM matrices. |
| **Part 20** | Prodigal · VIBRANT · eggNOG-mapper · AMRFinderPlus · DIAMOND | **Functional annotation:** ORF prediction; AMGs (**VIBRANT**); orthologs (**eggNOG-mapper** v2.1.12 / EggNOG 5.0); ARGs (**AMRFinderPlus**, **SARG** via DIAMOND); VFs (**VFDB** via DIAMOND). |
| **Part 21** | DIAMOND · MMseqs2 · `sbatch/diamond_ICTV.sbatch` · `sbatch/diamond_IMGVR.sbatch` | **Novelty assessment:** nucleotide and protein similarity to **ICTV** and **IMG/VR v4**; cross-study clustering against published porcine catalogues. |
| **Part 22** | `sbatch/vcontact3_gsvc.sbatch` | **Phylogenetic analysis (optional):** vContact3 network for selected viral families (Fig. 2b). |

**Output:** GSVC comprising 108,763 high-quality and 115,842 medium-quality vOTUs (224,605 total).

---

### Statistical analysis and visualization

Downstream ecological analyses, statistical tests, and figure generation are implemented in **R**, **Python**, **Shell**, and **Perl**. These scripts are placed in the **`Scripts`** directory, organized by main figure:

| Directory | Content |
|-----------|---------|
| `Scripts/Fig01_catalogue/` | Catalogue overview, novelty, read capture (Fig. 1) |
| `Scripts/Fig02_taxonomy/` | Taxonomic composition and phylogeny (Fig. 2) |
| `Scripts/Fig03_biogeography/` | Biogeography, PERMANOVA, distance-decay, neutral model (Fig. 3) |
| `Scripts/Fig04_host_phage/` | Host range, lifestyle, Kill-the-Winner dynamics (Fig. 4) |
| `Scripts/Fig05_core_enterotype/` | Prevalence, core virome, enterotypes (Fig. 5) |
| `Scripts/Fig06_function/` | ARG / AMG / functional coupling (Fig. 6) |
| `Scripts/Fig07_intervention/` | Post-weaning diarrhea cohort re-analysis (Fig. 7) |

All analysis-ready input tables for figure reproduction are in **`Pre-processed_Files/`** (see `Pre-processed_Files/MANIFEST.tsv` for a full inventory).

Large data files (Count matrix, vOTU FASTA, iPHoP predictions, CRISPR MAGs) are hosted on **Zenodo** ([10.5281/zenodo.20579313](https://doi.org/10.5281/zenodo.20579313)) — see **`large_data/README.md`**.

---

## REPOSITORY STRUCTURE

```text
GSVC_github/
├── Pipeline/              # Upstream metagenomic workflow
│   ├── modules/           # Shell modules
│   ├── sbatch/            # SLURM job templates
│   ├── script/            # Python/R helper scripts
│   └── pipeline.txt       # End-to-end command log
├── Pre-processed_Files/   # Metadata and analysis-ready tables
├── Scripts/                 # Figure reproduction and statistics
├── large_data/              # Zenodo download instructions
├── README.md
└── LICENSE                  # CC BY 4.0
```

---

## DATA AVAILABILITY

| Resource | Location |
|----------|----------|
| Scripts and pre-processed tables | This GitHub repository |
| vOTU sequences, Count matrix, iPHoP table, CRISPR MAGs | [Zenodo 10.5281/zenodo.20579313](https://doi.org/10.5281/zenodo.20579313) (see `large_data/README.md`) |
| Source metagenomes | NCBI SRA / ENA / CNGB (Table S1 in manuscript) |

---

## LICENSE

Released under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE).

## CONTACT

Questions and bug reports are welcome via GitHub Issues.
