#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="合并样本多次测序的fq路径为逗号分隔格式")
    parser.add_argument('-i', '--input', required=True, help='输入样本原始表格路径（含多行）')
    parser.add_argument('-o', '--output', required=True, help='输出整理后的样本表格路径')
    args = parser.parse_args()

    # 读取原始数据
    df = pd.read_csv(args.input, sep="\t")

    # 检查必须列
    for col in ['sample_id', 'bioproject', 'fq1', 'fq2']:
        if col not in df.columns:
            raise ValueError(f"缺少必要列: {col}")

    # 合并 fq1/fq2
    df_agg = df.groupby(["sample_id", "bioproject"]).agg({
        "fq1": lambda x: ",".join(sorted(x)),
        "fq2": lambda x: ",".join(sorted(x))
    }).reset_index()

    # 输出为 tsv
    df_agg.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] 整理完成，输出文件：{args.output}")

if __name__ == "__main__":
    main()
