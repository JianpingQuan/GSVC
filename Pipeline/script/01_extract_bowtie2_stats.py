import os
import csv
import argparse

def parse_bowtie2_log(log_file):
    stats = {
        'total_reads': None,
        'concordant_0': None,
        'concordant_1': None,
        'concordant_multi': None,
        'unmapped_mates': None,
        'single_align_1': None,
        'single_align_multi': None,
        # 'overall_alignment_rate' 去掉了
    }

    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            try:
                if "reads; of these:" in line and stats['total_reads'] is None:
                    stats['total_reads'] = int(line.split()[0])
                elif "aligned concordantly 0 times" in line and "pairs" not in line:
                    stats['concordant_0'] = int(line.split()[0])
                elif "aligned concordantly exactly 1 time" in line:
                    stats['concordant_1'] = int(line.split()[0])
                elif "aligned concordantly >1 times" in line:
                    stats['concordant_multi'] = int(line.split()[0])
                elif "aligned 0 times" in line and "mates make up the pairs" not in line:
                    stats['unmapped_mates'] = int(line.split()[0])
                elif "aligned exactly 1 time" in line and stats['single_align_1'] is None:
                    stats['single_align_1'] = int(line.split()[0])
                elif "aligned >1 times" in line and stats['single_align_multi'] is None:
                    stats['single_align_multi'] = int(line.split()[0])
            except (IndexError, ValueError):
                print(f"[警告] 跳过格式不标准的行: {line}")

    return stats


def extract_all_samples(nohost_dir, bioproject):
    output = []
    samples = [s for s in os.listdir(nohost_dir) if os.path.isdir(os.path.join(nohost_dir, s))]

    for sample in samples:
        log_path = os.path.join(nohost_dir, sample, "bowtie2.log")
        if not os.path.exists(log_path):
            print(f"[跳过] 未找到日志文件: {log_path}")
            continue
        if os.path.getsize(log_path) == 0:
            print(f"[跳过] 日志文件为空: {log_path}")
            continue

        try:
            stats = parse_bowtie2_log(log_path)
            stats["sample"] = sample
            output.append(stats)
        except Exception as e:
            print(f"[错误] 解析 {log_path} 失败: {e}")

    return output


def write_csv(stats_list, bioproject, output_dir):
    output_file = os.path.join(output_dir, f"{bioproject}_bowtie2.csv")
    fieldnames = [
        "sample",
        "total_reads",
        "concordant_0",
        "concordant_1",
        "concordant_multi",
        "unmapped_mates",
        "single_align_1",
        "single_align_multi",
    ]

    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in stats_list:
            writer.writerow(row)

    print(f"[完成] 写入统计结果到：{output_file}")


def main():
    parser = argparse.ArgumentParser(description="Extract bowtie2 log statistics")
    parser.add_argument("--bioproject", required=True, help="Bioproject ID, e.g. PRJNA735412")
    parser.add_argument("--nohost_dir", required=True, help="Path to bowtie2 log root folder")
    args = parser.parse_args()

    stats = extract_all_samples(args.nohost_dir, args.bioproject)
    write_csv(stats, args.bioproject, args.nohost_dir)


if __name__ == "__main__":
    main()

