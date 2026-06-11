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
mapped_prefix="${outdir}/${sample}_map"
logfile="${outdir}/bowtie2.log"

# 运行bowtie2并同时输出日志
bowtie2 -x $GENOME_INDEX -1 $fq1 -2 $fq2 \
  --un-conc-gz ${unmapped_prefix}_%.fq.gz \
  --al-conc-gz ${mapped_prefix}_%.fq.gz \
  -S /dev/null -p $THREADS \
  2> $logfile

# 提取比对率并保存
rate=$(grep "overall alignment rate" $logfile | awk '{print $(NF-1) $NF}')  # e.g., "95.24%"
echo -e "${sample}\t${bioproject}\t${rate}" >> "${OUTDIR}/mapping_rate.tsv"
