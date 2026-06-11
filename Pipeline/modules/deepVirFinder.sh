#!/bin/bash

sample=$1
bioproject=$2


# 路径设置
OUTDIR="/work/home/acq6mmxy68/Public/results"
megahit_outdir="${OUTDIR}/${bioproject}/megahit/${sample}"
contig="${megahit_outdir}/${sample}.contigs.fa"
dvf_outdir="${OUTDIR}/${bioproject}/dvf/${sample}"

# 创建输出目录
mkdir -p "$dvf_outdir"

# 执行deepVirFinder
python /work/home/acq6mmxy68/software/DeepVirFinder/dvf.py -i "$contig" -o "${dvf_outdir}" -l 10000 -c 32
