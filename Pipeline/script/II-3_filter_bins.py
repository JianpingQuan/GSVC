import ast
import os
import sys

if len(sys.argv) != 4:
    print("Usage: python filter_bins.py <bin_stats_file> <bins_dir> <output_dir>")
    sys.exit(1)

bin_stats_file = sys.argv[1]
bins_dir = sys.argv[2]
output_dir = sys.argv[3]

os.makedirs(output_dir, exist_ok=True)

# 新增一个统计文件路径
summary_file = os.path.join(output_dir, "bins_completeness_contamination.tsv")

with open(bin_stats_file, 'r') as f, open(summary_file, 'w') as summary_out:
    # 写入表头
    summary_out.write("Bin_ID\tCompleteness\tContamination\n")

    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            bin_id, info_str = line.split('\t', 1)
            info = ast.literal_eval(info_str)
            comp = info.get('Completeness', 0)
            cont = info.get('Contamination', 100)

            # 写入所有bin的统计信息
            summary_out.write(f"{bin_id}\t{comp}\t{cont}\n")

            # 过滤并复制符合条件的bins
            if comp >= 90 and cont <= 5:
                src = os.path.join(bins_dir, f"{bin_id}.fa")
                dst = os.path.join(output_dir, f"{bin_id}.fa")
                if os.path.exists(src):
                    os.system(f"cp {src} {dst}")
                else:
                    print(f"[Warning] File not found: {src}")
        except Exception as e:
            print(f"[Error] Failed to process line: {line}")
            print(e)
