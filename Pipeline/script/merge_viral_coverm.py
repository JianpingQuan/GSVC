import os
import sys
import pandas as pd

def merge_viral_data(bioproject_id):
    # 基础路径（根据你提供的最新路径）
    base_path = f'/work/home/acq6mmxy68/Public/results/{bioproject_id}/viral_coverM'
    out_tpm = os.path.join(base_path, f'merged_{bioproject_id}_TPM.tsv')
    out_count = os.path.join(base_path, f'merged_{bioproject_id}_Count.tsv')

    if not os.path.exists(base_path):
        print(f"错误: 找不到路径 {base_path}")
        return

    tpm_dfs = []
    count_dfs = []
    
    # 获取该目录下所有子文件夹并排序
    subdirs = sorted([d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))])

    print(f"开始处理项目: {bioproject_id}")
    print(f"搜索路径: {base_path}")

    for subdir in subdirs:
        # 寻找 output_count.tsv 文件
        file_path = os.path.join(base_path, subdir, 'output_count.tsv')
        
        if os.path.exists(file_path):
            try:
                # 读取文件，第一列(Contig)作为索引
                # sep='\t' 处理制表符，index_col=0 将第一列设为行名
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                
                # 提取 Read Count (文件中的第二列，即 iloc 索引 0)
                count_series = df.iloc[:, 0].copy()
                count_series.name = f"{subdir}_Count"
                count_dfs.append(count_series)
                
                # 提取 TPM (文件中的第三列，即 iloc 索引 1)
                tpm_series = df.iloc[:, 1].copy()
                tpm_series.name = f"{subdir}_TPM"
                tpm_dfs.append(tpm_series)
                
                print(f"  [成功] 已加载: {subdir}")
            except Exception as e:
                print(f"  [跳过] {subdir} 出现错误: {e}")
        else:
            # 如果某个子文件夹里没有这个文件，打印提示
            print(f"  [跳过] {subdir} 目录下未找到 output_count.tsv")

    # 执行合并
    if tpm_dfs and count_dfs:
        print("\n正在执行矩阵合并与空值填充...")
        
        # 合并 Count 并保存
        final_count = pd.concat(count_dfs, axis=1, join='outer').fillna(0)
        final_count.to_csv(out_count, sep='\t')
        print(f">>> Count 矩阵已保存至: {out_count}")
        
        # 合并 TPM 并保存
        final_tpm = pd.concat(tpm_dfs, axis=1, join='outer').fillna(0)
        final_tpm.to_csv(out_tpm, sep='\t')
        print(f">>> TPM 矩阵已保存至: {out_tpm}")
    else:
        print("未发现任何有效数据文件，请检查路径。")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使用方法: python3 merge_viral_coverm.py <BIOPROJECT_ID>")
        print("示例: python3 merge_viral_coverm.py CNP0005498")
    else:
        merge_viral_data(sys.argv[1])
