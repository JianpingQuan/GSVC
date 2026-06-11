#!/bin/bash

# 输入：MEGAHIT 总输出目录
megahit_dir=$1

# 检查 QUAST 是否安装
if ! command -v quast.py &> /dev/null; then
  echo "[ERROR] quast.py not found. Please activate your conda or module environment."
  exit 1
fi

# 遍历每个样本子目录
for sample_dir in "$megahit_dir"/*/; do
  sample=$(basename "$sample_dir")
  contig1="${sample_dir}/final.contigs.fa"
  contig2="${sample_dir}/${sample}.contigs.fa"

  if [ -f "$contig1" ]; then
    contig_file="$contig1"
  elif [ -f "$contig2" ]; then
    contig_file="$contig2"
  else
    echo "[WARNING] No contig file found for $sample"
    continue
  fi

  outdir="${sample_dir}/quast_report"
  mkdir -p "$outdir"

  echo "[INFO] Running QUAST for $sample..."
  quast.py "$contig_file" -o "$outdir" -t 8 --min-contig 1000 > "${outdir}/quast.log" 2>&1
done
