#!/bin/bash
set -e

# 接收外部输入
bioproject="$1"
sample="$2"

# 检查参数
if [ -z "$bioproject" ] || [ -z "$sample" ]; then
  echo "[Usage: bash rename_bins.sh <bioproject> <sample>"
  exit 1
fi

# 输入与输出路径
refined_dir="/work/home/acq6mmxy68/Public/results/${bioproject}/checkm/${sample}/refined_bins"
output_dir="/work/home/acq6mmxy68/Public/results/${bioproject}/drep/drep_input"
mkdir -p "$output_dir"

# 遍历 refined_bins 中的 bin 文件
for bin_file in "$refined_dir"/bin.*.fa; do
  [ -e "$bin_file" ] || continue

  bin_id=$(basename "$bin_file" .fa)  # bin.104
  new_name="${bioproject}_${sample}.${bin_id}.fa"
  output_file="${output_dir}/${new_name}"

  # 替换 header 并输出
  awk -v prefix="${bioproject}_${sample}_" '/^>/ {sub(/^>/, ">"prefix); print; next} {print}' "$bin_file" > "$output_file"

  echo "[✔] Renamed: $(basename "$bin_file") → $(basename "$output_file")"
done
