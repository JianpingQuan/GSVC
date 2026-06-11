import os
import sys
import pandas as pd

def merge_data(target_dir):
    # 获取项目 ID（取路径的最后一级作为文件名标识）
    project_id = os.path.basename(target_dir.rstrip('/'))
    
    out_tpm = os.path.join(target_dir, f'merged_{project_id}_TPM.tsv')
    out_count = os.path.join(target_dir, f'merged_{project_id}_Count.tsv')

    if not os.path.exists(target_dir):
        print(f"错误: 找不到路径 {target_dir}")
        return

    tpm_dfs = []
    count_dfs = []
    
    # 获取子文件夹
    subdirs = sorted([d for d in os.listdir(target_dir) if os.path.isdir(os.path.join(target_dir, d))])

    print(f"开始处理路径: {target_dir}")

    for subdir in subdirs:
        file_path = os.path.join(target_dir, subdir, 'output_count.tsv')
        
        if os.path.exists(file_path):
            try:
                # 读取文件，第一列(Contig)为索引
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                
                # 提取 Read Count (第二列)
                count_series = df.iloc[:, 0].copy()
                count_series.name = f"{subdir}_Count"
                count_dfs.append(count_series)
                
                # 提取 TPM (第三列)
                tpm_series = df.iloc[:, 1].copy()
                tpm_series.name = f"{subdir}_TPM"
                tpm_dfs.append(tpm_series)
                
                print(f"  [OK] {subdir}")
            except Exception as e:
                print(f"  [ERROR] {subdir}: {e}")

    # 合并并保存
    if tpm_dfs and count_dfs:
        print("\n正在生成大矩阵...")
        
        # 合并 Count
        pd.concat(count_dfs, axis=1, join='outer').fillna(0).to_csv(out_count, sep='\t')
        print(f">>> Count 矩阵已保存: {out_count}")
        
        # 合并 TPM
        pd.concat(tpm_dfs, axis=1, join='outer').fillna(0).to_csv(out_tpm, sep='\t')
        print(f">>> TPM 矩阵已保存: {out_tpm}")
    else:
        print("未发现 output_count.tsv 文件。")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使用方法: python3 universal_merge_coverm.py <TARGET_DIRECTORY_PATH>")
    else:
        merge_data(sys.argv[1])
