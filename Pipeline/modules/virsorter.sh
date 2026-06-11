#!/bin/bash

sample=$1
bioproject=$2


# 路径设置
OUTDIR="/work/home/acq6mmxy68/Public/results"
megahit_outdir="${OUTDIR}/${bioproject}/megahit/${sample}"
contig="${megahit_outdir}/${sample}.contigs.fa"
virsorter_outdir="${OUTDIR}/${bioproject}/virsorter/${sample}"

# 创建输出目录
mkdir -p "$virsorter_outdir"

# 执行virsorter2
virsorter run -w "${virsorter_outdir}" -i "$contig" --min-length 10000 --include-groups dsDNAphage,NCLDV,ssDNA --high-confidence-only -j 32 all
