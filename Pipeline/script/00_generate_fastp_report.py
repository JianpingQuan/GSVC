import os
import json
import csv
import sys

def parse_fastp_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    summary = data.get('summary', {})
    before = summary.get('before_filtering', {})
    after = summary.get('after_filtering', {})

    return {
        'total_reads_before': before.get('total_reads'),
        'total_reads_after': after.get('total_reads'),
        'total_bases_before': before.get('total_bases'),
        'total_bases_after': after.get('total_bases'),
        'q20_rate_before': before.get('q20_rate'),
        'q20_rate_after': after.get('q20_rate'),
        'q30_rate_before': before.get('q30_rate'),
        'q30_rate_after': after.get('q30_rate'),
        'gc_content_before': before.get('gc_content'),
        'gc_content_after': after.get('gc_content'),
    }

def generate_summary_from_clean_dir(bioproject, clean_data_dir):
    results = []
    sample_dirs = [d for d in os.listdir(clean_data_dir) if os.path.isdir(os.path.join(clean_data_dir, d))]

    for sample in sample_dirs:
        sample_path = os.path.join(clean_data_dir, sample)
        json_files = [f for f in os.listdir(sample_path) if f.endswith('.json')]

        if not json_files:
            print(f"[警告] {sample_path} 中未找到 JSON 文件")
            continue

        json_file_path = os.path.join(sample_path, json_files[0])
        try:
            fastp_data = parse_fastp_json(json_file_path)
            fastp_data['sample'] = sample
            fastp_data['bioproject'] = bioproject
            results.append(fastp_data)
        except Exception as e:
            print(f"[错误] 解析 {json_file_path} 失败: {e}")

    if results:
        output_file = os.path.join(clean_data_dir, f"{bioproject}_fastp_summary.csv")
        fieldnames = ['sample', 'bioproject'] + [k for k in results[0] if k not in ('sample', 'bioproject')]

        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"[完成] 汇总表已生成: {output_file}")
    else:
        print("[提示] 没有成功解析任何样本 JSON 文件。")

# 如果作为命令行工具使用
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python generate_fastp_report.py <bioproject_id> <clean_data_dir>")
        print("示例: python generate_fastp_report.py PRJNA735412 /path/to/clean")
        sys.exit(1)

    bioproject_id = sys.argv[1]
    clean_data_dir = sys.argv[2]

    generate_summary_from_clean_dir(bioproject_id, clean_data_dir)
