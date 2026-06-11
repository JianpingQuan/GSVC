#!/bin/bash
set -euo pipefail

# 输入参数：MEGAHIT输出路径
megahit_dir=$1

# 检查路径存在
if [ ! -d "$megahit_dir" ]; then
  echo "[ERROR] Directory not found: $megahit_dir"
  exit 1
fi

# 输出统计结果文件
summary_file="${megahit_dir}/contig_stats_summary.tsv"
echo -e "Sample\tContig_Count\tTotal_Length\tMax_Length\tAvg_Length\tN50" > "$summary_file"

# 遍历子目录（假设子目录名为样本名）
for sample_dir in "$megahit_dir"/*; do
  if [ -d "$sample_dir" ]; then
    sample=$(basename "$sample_dir")
    contig_file="${sample_dir}/${sample}.contigs.fa"

    if [ ! -f "$contig_file" ]; then
      echo "[WARN] Skipped missing contig file: $contig_file"
      continue
    fi

    awk -v sample="$sample" '
      BEGIN {
        seq = ""; total_len = 0; max_len = 0; contig_count = 0;
      }
      /^>/ {
        if (length(seq) > 0) {
          len = length(seq);
          total_len += len;
          contig_lens[contig_count++] = len;
          if (len > max_len) max_len = len;
        }
        seq = "";
        next;
      }
      {
        seq = seq $0;
      }
      END {
        if (length(seq) > 0) {
          len = length(seq);
          total_len += len;
          contig_lens[contig_count++] = len;
          if (len > max_len) max_len = len;
        }

        avg_len = (contig_count > 0) ? total_len / contig_count : 0;

        asort(contig_lens, sorted_lens, "@val_num_desc");
        sum = 0; n50 = 0;
        for (i = 1; i <= contig_count; i++) {
          sum += sorted_lens[i];
          if (sum >= total_len / 2) {
            n50 = sorted_lens[i];
            break;
          }
        }

        printf "%s\t%d\t%d\t%d\t%.2f\t%d\n", sample, contig_count, total_len, max_len, avg_len, n50;
      }
    ' "$contig_file" >> "$summary_file"
  fi
done

echo "[INFO] Summary saved to: $summary_file"
