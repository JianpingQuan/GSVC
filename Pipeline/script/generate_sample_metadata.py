#!/usr/bin/env python3
import os
import pandas as pd
import re
import argparse
from collections import defaultdict

def find_paired_files(data_dir):
    """递归查找并配对测序文件"""
    pattern_pairs = [
        # 匹配 _1.fq.gz 或 _2.fq.gz 模式
        (r'^(.+)_([12])\.fq\.gz$', 1),
        # 匹配 _1.fastq.gz 或 _2.fastq.gz 模式
        (r'^(.+)_([12])\.fastq\.gz$', 1),
        # 匹配 _f1.fq.gz 或 _r2.fq.gz 模式
        (r'^(.+)_[fr]([12])\.fq\.gz$', 1),
        # 匹配 .1.fq.gz 或 .2.fq.gz 模式
        (r'^(.+)\.([12])\.fq\.gz$', 1),
    ]
    
    file_groups = defaultdict(dict)
    
    # 递归遍历所有子目录
    for root, dirs, files in os.walk(data_dir):
        for filename in files:
            if not (filename.endswith('.fq.gz') or filename.endswith('.fastq.gz')):
                continue
                
            for pattern, ext_group in pattern_pairs:
                match = re.match(pattern, filename)
                if match:
                    prefix = match.group(1)
                    pair_num = match.group(2)
                    sample_id = os.path.basename(root)  # 使用子目录名作为样本ID
                    
                    file_path = os.path.join(root, filename)
                    file_groups[(sample_id, prefix)][pair_num] = file_path
                    break
    
    return file_groups

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Generate sample metadata TSV for a given BioProject")
    parser.add_argument("--bioproject", required=True, help="BioProject number (e.g. PRJNA735412)")
    parser.add_argument("--data_dir", required=True, help="Directory containing sample subdirectories")
    args = parser.parse_args()
    
    file_groups = find_paired_files(args.data_dir)
    
    records = []
    for (sample_id, prefix), files in file_groups.items():
        fq1 = files.get('1')
        fq2 = files.get('2')
        
        if fq1 and fq2:
            records.append({
                "sample_id": sample_id,
                "bioproject": args.bioproject,
                "fq1": fq1,
                "fq2": fq2
            })
        elif fq1:
            print(f"[WARN] Only fq1 found for {sample_id}: {fq1}")
        elif fq2:
            print(f"[WARN] Only fq2 found for {sample_id}: {fq2}")
    
    # 输出路径
    output_dir = "/work/home/acq6mmxy68/pipeline/config/sampleList"
    output_path = f"{output_dir}/sample_metadata_{args.bioproject}.tsv"
    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.DataFrame(records)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"[INFO] Metadata file written to: {output_path}")
    print(f"[INFO] Found {len(records)} valid sample pairs")

if __name__ == "__main__":
    main()
