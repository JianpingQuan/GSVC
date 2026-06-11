#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from collections import defaultdict
import re

def find_clean_paired_files(data_dir):
    """递归查找并配对 *_1.clean.fq.gz 和 *_2.clean.fq.gz 文件"""
    file_groups = defaultdict(dict)
    pattern = r'^(.+)_([12]).clean\.fq\.gz$'

    for root, dirs, files in os.walk(data_dir):
        for filename in files:
            if not filename.endswith('.fq.gz'):
                continue

            match = re.match(pattern, filename)
            if not match:
                continue

            prefix = match.group(1)
            pair_num = match.group(2)
            sample_id = os.path.basename(root)
            file_path = os.path.join(root, filename)
            file_groups[(sample_id, prefix)][pair_num] = file_path

    return file_groups

def main():
    parser = argparse.ArgumentParser(description="Generate metadata for clean fq.gz pairs")
    parser.add_argument("--bioproject", required=True, help="BioProject ID (e.g. PRJNAxxxxxx)")
    parser.add_argument("--data_dir", required=True, help="Root directory containing sample subfolders")
    args = parser.parse_args()

    file_groups = find_clean_paired_files(args.data_dir)
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
        else:
            print(f"[WARN] Incomplete pair for {sample_id}/{prefix}: fq1={bool(fq1)} fq2={bool(fq2)}")

    output_dir = "/work/home/acq6mmxy68/pipeline/config/sampleList/"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"sample_metadata_{args.bioproject}_clean.tsv")

    df = pd.DataFrame(records)
    df.to_csv(output_path, sep="\t", index=False)

    print(f"[INFO] Metadata file written to: {output_path}")
    print(f"[INFO] Found {len(records)} valid clean sample pairs")

if __name__ == "__main__":
    main()
