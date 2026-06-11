# Large data files (Zenodo)

Files too large for GitHub are hosted on Zenodo:

**DOI:** [10.5281/zenodo.20579313](https://doi.org/10.5281/zenodo.20579313)  
**Record:** https://zenodo.org/records/20579313

Download all required archives, verify MD5 checksums, extract, and place files as described below.

---

## File inventory

| File | Size | MD5 | Description |
|------|------|-----|-------------|
| `FINAL_all_projects_Count_matrix.zip` | 380.03 MB | `0499d2ee519bff17daa13f5e15c13c6f` | vOTU × sample raw read counts (CoverM); unzip to `.tsv` |
| `high_quality_viral.Final.fasta.gz` | 1.76 GB | `f35e3aa440779618129af3278cf2cff0` | High-quality vOTU representative genomes (≥90% complete) |
| `medium_quality_viral.Final.fasta.gz` | 1.28 GB | `80c4b70264427c1853e05fd49f6f5576` | Medium-quality vOTU representative genomes (50–90% complete) |
| `high_quality_Host_prediction_to_genome_m90.csv` | 62.67 MB | `b74c82009871ee8d01b2ea42b2ea5b0e` | iPHoP host predictions (genome level, medium+ confidence) |
| `CRISPR_MAGs_01.tar.gz` | 956.02 MB | `ba380b5380ce714c375de5a61df491af` | Bacterial MAGs for CRISPR spacer calling (part 1/10) |
| `CRISPR_MAGs_02.tar.gz` | 952.43 MB | `ca940f1c8025a99c1d4cd1c804a63cf4` | Bacterial MAGs for CRISPR spacer calling (part 2/10) |
| `CRISPR_MAGs_03.tar.gz` | 947.08 MB | `4dd539a006273ece4208d95e4edffb10` | Bacterial MAGs for CRISPR spacer calling (part 3/10) |
| `CRISPR_MAGs_04.tar.gz` | 949.24 MB | `9d70b559befa41f682c4643eb910ce97` | Bacterial MAGs for CRISPR spacer calling (part 4/10) |
| `CRISPR_MAGs_05.tar.gz` | 966.72 MB | `b74b3016fbf88a0c4730ada0512a5ca5` | Bacterial MAGs for CRISPR spacer calling (part 5/10) |
| `CRISPR_MAGs_06.tar.gz` | 944.40 MB | `f4f177d46e44ae5d4af371704b53d86b` | Bacterial MAGs for CRISPR spacer calling (part 6/10) |
| `CRISPR_MAGs_07.tar.gz` | 942.62 MB | `e786cbf86d3ab85f60e16c4086def4a2` | Bacterial MAGs for CRISPR spacer calling (part 7/10) |
| `CRISPR_MAGs_08.tar.gz` | 945.29 MB | `06347bad858b79b1fd4a7563221298ee` | Bacterial MAGs for CRISPR spacer calling (part 8/10) |
| `CRISPR_MAGs_09.tar.gz` | 946.10 MB | `13887aa5397041d1c1d7a426eac068fe` | Bacterial MAGs for CRISPR spacer calling (part 9/10) |
| `CRISPR_MAGs_10.tar.gz` | 294.6  MB | `e1b4bceed5e503102ffc7bfbd6b26c0b` | Bacterial MAGs for CRISPR spacer calling (part 10/10) |

---

## Download and setup

```bash
# Example: download via Zenodo record page or wget
# https://zenodo.org/records/20579313

# Verify checksum (Linux / macOS)
md5sum FINAL_all_projects_Count_matrix.zip
# expected: 0499d2ee519bff17daa13f5e15c13c6f

# Extract Count matrix to repo root
unzip FINAL_all_projects_Count_matrix.zip -d ../
# or: unzip -p FINAL_all_projects_Count_matrix.zip > ../FINAL_all_projects_Count_matrix.tsv

# Extract vOTU FASTA
gunzip -k high_quality_viral.Final.fasta.gz medium_quality_viral.Final.fasta.gz

# Extract CRISPR MAG archives
mkdir -p ../Database/CRISPR_MAGs
for f in CRISPR_MAGs_*.tar.gz; do tar -xzf "$f" -C ../Database/CRISPR_MAGs/; done
```

```powershell
# Windows PowerShell — verify MD5
Get-FileHash FINAL_all_projects_Count_matrix.zip -Algorithm MD5

# Extract zip to repo root
Expand-Archive FINAL_all_projects_Count_matrix.zip -DestinationPath ..
```

These large files are listed in `.gitignore` and must not be committed to Git.

---

## Required for which analyses?

| File | Required by |
|------|-------------|
| `FINAL_all_projects_Count_matrix.tsv` | Fig. 3 (biogeography), Fig. 5 (prevalence, core, enterotype), PERMANOVA |
| `high_quality_viral.Final.fasta` + `medium_quality_viral.Final.fasta` | Custom re-profiling of new metagenomes; CRISPR spacer matching |
| `high_quality_Host_prediction_to_genome_m90.csv` | Fig. 4 host–phage analyses (full iPHoP table) |
| `CRISPR_MAGs_01–10.tar.gz` | Re-running CRISPR spacer identification (`Pipeline/script/identify_crispr.py`) |

---

## Re-profiling new samples

To quantify GSVC vOTUs in new pig gut metagenomes:

1. Quality-filter and assemble reads (see `Pipeline/` and manuscript Methods).  
2. Concatenate or index `high_quality_viral.Final.fasta` and `medium_quality_viral.Final.fasta` as the reference database.  
3. Map reads with CoverM (or equivalent; ≥95% nucleotide identity recommended for strict matching).  
4. Use the same detection rule as the manuscript: **count > 0** defines presence; relative abundance = vOTU count / total vOTU counts per sample.

---