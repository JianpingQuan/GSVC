#!/bin/bash

sample=$1
bioproject=$2

paired_path="/work/home/acq6mmxy68/Public/results/${bioproject}/nohost/${sample}"
genome_path="/work/home/acq6mmxy68/"
out_path="/work/home/acq6mmxy68/Public/results/${bioproject}/viral_coverM/${sample}"

mkdir -p $out_path

coverm contig \
  --coupled "${paired_path}/${sample}_unmap_1.fq.gz" "${paired_path}/${sample}_unmap_2.fq.gz" \
  --reference ${genome_path}/PGV.fna \
  --min-read-aligned-percent 75 \
  --min-read-percent-identity 95 \
  --min-covered-fraction 75 \
  -t 32 \
  --methods tpm \
  -o $out_path/output_coverm_mapped_study1.tsv


#coverm contig \
#  --coupled "${paired_path}/${sample}_unmap_1.fq.gz" "${paired_path}/${sample}_unmap_2.fq.gz" \
#  --reference ${genome_path}/high_medium_quality_viral.Final.fasta \
#  --min-read-aligned-percent 75 \
#  --min-read-percent-identity 95 \
#  -t 32 \
#  --methods count tpm\
 # -o $out_path/output_count.tsv
