#!/bin/bash

source /work/home/acq6mmxy68/pipeline/config/config.param

sample=$1
bioproject=$2
fq1=$3   # 支持逗号分隔多个 fq.gz
fq2=$4

outdir="${OUTDIR}/${bioproject}/clean/${sample}"
mkdir -p "$outdir"

# 临时合并后的 fq 文件路径
fq1_merged="${outdir}/${sample}_merged_1.fq.gz"
fq2_merged="${outdir}/${sample}_merged_2.fq.gz"

# 判断是否是多个文件，如果是就合并
if [[ "$fq1" == *,* ]]; then
  echo "[INFO] Detected multiple R1 files. Merging..."
  cat $(echo $fq1 | tr ',' ' ') > "$fq1_merged"
else
  cp "$fq1" "$fq1_merged"
fi

if [[ "$fq2" == *,* ]]; then
  echo "[INFO] Detected multiple R2 files. Merging..."
  cat $(echo $fq2 | tr ',' ' ') > "$fq2_merged"
else
  cp "$fq2" "$fq2_merged"
fi

# fastp QC
fastp -i "$fq1_merged" -I "$fq2_merged" \
  -o "${outdir}/${sample}_1.clean.fq.gz" \
  -O "${outdir}/${sample}_2.clean.fq.gz" \
  -w $THREADS \
  -h "${outdir}/fastp.html" \
  -j "${outdir}/fastp.json"
