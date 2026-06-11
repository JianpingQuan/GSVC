import pandas as pd
import glob
import os

# 1. 设置搜索路径
search_path = "/work/home/acq6mmxy68/Public/results/*/coverM/*/output_coverm.tsv"
files = sorted(glob.glob(search_path))

print(f"找到 {len(files)} 个文件，开始集成...")

data_dict = {}

for f_path in files:
    # 提取 Sample ID (取文件名所在目录的名字)
    sample_id = f_path.split('/')[-2]
    
    try:
        # 读取数据，设置第一列 Genome 为索引
        df = pd.read_csv(f_path, sep='\t', index_col=0)
        
        # --- 关键修改：剔除 unmapped 行 ---
        # 使用 errors='ignore' 确保即使某文件没有 unmapped 行也不会报错
        df = df.drop('unmapped', errors='ignore')
        
        # 寻找包含 "TPM" 字样的列名（解决表头带样本名导致的匹配问题）
        tpm_col = [c for c in df.columns if 'TPM' in c]
        
        if not tpm_col:
            print(f"警告: 文件 {sample_id} 中未找到 TPM 列")
            continue
            
        # 提取 TPM 列并转换为数值型，强制将无法转换的(如NA)变为 NaN
        data_dict[sample_id] = pd.to_numeric(df[tpm_col[0]], errors='coerce')
        
    except Exception as e:
        print(f"处理样本 {sample_id} 时发生错误: {e}")

# 2. 合并所有样本数据
# 使用 outer join 确保所有基因组 ID 都能保留，缺失值自动补 NaN
combined_matrix = pd.concat(data_dict.values(), axis=1, keys=data_dict.keys(), sort=True)

# 3. 后处理：将所有的 NaN (缺失或之前的 NA) 填充为 0.0
combined_matrix.fillna(0.0, inplace=True)
combined_matrix.index.name = 'Genome'

# 4. 保存结果
output_file = "merged_tpm_clean.tsv"
combined_matrix.to_csv(output_file, sep='\t')

print("-" * 30)
print(f"合并成功！")
print(f"最终矩阵规模: {combined_matrix.shape[0]} 行 (Genomes) x {combined_matrix.shape[1]} 列 (Samples)")
print(f"文件已保存至: {output_file}")
