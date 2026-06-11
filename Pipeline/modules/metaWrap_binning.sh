#!/bin/bash

sample=$1
bioproject=$2

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1


#
output_dir="/work/home/acq6mmxy68/Public/results/${bioproject}/binningWrap/${sample}/bin"
assembly="/work/home/acq6mmxy68/Public/results/${bioproject}/megahit/${sample}/${sample}.contigs.fa"
read1="/work/home/acq6mmxy68/Public/results/${bioproject}/nohost/${sample}/${sample}_unmap_1.fastq"
read2="/work/home/acq6mmxy68/Public/results/${bioproject}/nohost/${sample}/${sample}_unmap_2.fastq"


# 创建输出目录
mkdir -p "$output_dir"

# 执行 metaWrap 进行 binning
metawrap binning \
  -o "$output_dir" \
  -t 32 -m 4 \
  -a "$assembly" \
  --metabat2 --maxbin2 \
  "$read1" "$read2"
