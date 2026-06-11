#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
过滤高质量和中等质量的MAG，并分别保存到不同文件夹。
高质量: completeness >= 90, contamination <= 5
中等质量: 50 < completeness < 90, 5 < contamination < 10
"""

import ast
import os
import sys
import shutil

if len(sys.argv) != 4:
    print("Usage: python filter_bins.py <bin_stats_file> <bins_dir> <output_dir>")
    sys.exit(1)

bin_stats_file = sys.argv[1]
bins_dir = sys.argv[2]
output_dir = sys.argv[3]

# === 创建输出目录结构 ===
high_dir = os.path.join(output_dir, "high_quality_bins")
medium_dir = os.path.join(output_dir, "medium_quality_bins")
os.makedirs(high_dir, exist_ok=True)
os.makedirs(medium_dir, exist_ok=True)

summary_file = os.path.join(output_dir, "bins_completeness_contamination.tsv")

count_high = 0
count_medium = 0
count_low = 0

with open(bin_stats_file, 'r') as f, open(summary_file, 'w') as summary_out:
    summary_out.write("Bin_ID\tCompleteness\tContamination\tQuality\n")

    for line in f:
        line = line.strip()
        if not line:
            continue

        try:
            bin_id, info_str = line.split('\t', 1)
            info = ast.literal_eval(info_str)
            comp = float(info.get('Completeness', 0))
            cont = float(info.get('Contamination', 100))
            quality = "low"

            src = os.path.join(bins_dir, f"{bin_id}.fa")
            dst = None

            # 判断质量等级
            if comp >= 90 and cont <= 5:
                quality = "high"
                dst = os.path.join(high_dir, f"{bin_id}.fa")
                count_high += 1
            elif 50 < comp < 90 and 5 < cont < 10:
                quality = "medium"
                dst = os.path.join(medium_dir, f"{bin_id}.fa")
                count_medium += 1
            else:
                count_low += 1

            # 写入统计表
            summary_out.write(f"{bin_id}\t{comp}\t{cont}\t{quality}\n")

            # 拷贝文件
            if dst:
                if os.path.exists(src):
                    shutil.copy(src, dst)
                else:
                    print(f"[Warning] File not found: {src}")

        except Exception as e:
            print(f"[Error] Failed to process line: {line}")
            print(e)

# === 汇总输出 ===
total = count_high + count_medium + count_low
print("\n=== Filter Summary ===")
print(f"Total bins processed: {total}")
print(f"High-quality MAGs (>=90% completeness, <=5% contamination): {count_high}")
print(f"Medium-quality MAGs (50-90% completeness, 5-10% contamination): {count_medium}")
print(f"Low-quality or excluded: {count_low}")
print(f"\nResults saved in: {output_dir}\n")
