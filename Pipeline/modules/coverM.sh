#!/bin/bash

sample=$1
bioproject=$2

paired_path="/work/home/acq6mmxy68/Public/results/${bioproject}/nohost/${sample}"
genome_path="/work/home/acq6mmxy68/Public/results/00_CombinedAnalysis/allCleanMAGs_ANI95"
out_path="/work/home/acq6mmxy68/Public/results/${bioproject}/coverM/${sample}"

mkdir -p $out_path

coverm genome \
  --read1 "${paired_path}/${sample}_unmap_1.fq.gz" \
  --read2 "${paired_path}/${sample}_unmap_2.fq.gz" \
  --genome-fasta-files ${genome_path}/*.fa \
  -t 36 \
  -m mean relative_abundance covered_fraction tpm \
  -o $out_path/output_coverm.tsv
