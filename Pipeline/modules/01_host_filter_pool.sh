#!/bin/bash

set -e
source /work/home/acq6mmxy68/pipeline/config/config.param

sample=$1
bioproject=$2

fq1="${OUTDIR}/${bioproject}/clean/${sample}/${sample}_1.clean.fq.gz"
fq2="${OUTDIR}/${bioproject}/clean/${sample}/${sample}_2.clean.fq.gz"
outdir="${OUTDIR}/${bioproject}/nohost/${sample}"
mkdir -p $outdir

# 输出文件前缀
unmapped_prefix="${outdir}/${sample}_unmap"
logfile="${outdir}/bowtie2.log"

# 运行bowtie2并同时输出日志
bowtie2 -x $GENOME_INDEX -1 $fq1 -2 $fq2 \
  --un-conc-gz ${unmapped_prefix}_%.fq.gz \
  -S /dev/null -p $THREADS \
  2> $logfile
