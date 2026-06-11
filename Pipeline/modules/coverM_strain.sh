#!/bin/bash

sample=$1
bioproject=$2

paired_path="/work/home/acq6mmxy68/Public/results/${bioproject}/nohost/${sample}"
out_path="/work/home/acq6mmxy68/Public/results/${bioproject}/coverM/${sample}"

mkdir -p $out_path

coverm genome \
  --read1 "${paired_path}/${sample}_unmap_1.fq.gz" \
  --read2 "${paired_path}/${sample}_unmap_2.fq.gz" \
  -r /work/home/acq6mmxy68/dRep_strain_final_genomes/Minimap2_Index/all_MAGs_combined.mmi \
  --minimap2-reference-is-index \
  --genome-definition /work/home/acq6mmxy68/dRep_strain_final_genomes/genome_definition.tsv \
  -t 36 \
  --methods relative_abundance tpm \
  --min-read-aligned-percent 75 \
  --min-covered-fraction 0.1 \
  --min-read-percent-identity 99 \
  -o $out_path/output_coverm_strain.tsv
