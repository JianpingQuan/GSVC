import pandas as pd
import os

# 输出文件
output_file = "coverm_merged.tsv"
error_log = "error_files.txt"

df_final = None
error_files = []
success_count = 0

print("开始处理CoverM文件...")

with open("coverm.path.txt") as f:
    total_files = sum(1 for _ in f)
    f.seek(0)  # 回到文件开头
    
    for i, path in enumerate(f, 1):
        path = path.strip()
        if not path:  # 跳过空行
            continue

        # 显示进度
        if i % 100 == 0:
            print(f"处理进度: {i}/{total_files}")

        try:
            # 检查文件是否存在
            if not os.path.exists(path):
                error_files.append(f"{path} - 文件不存在")
                continue

            # 从路径提取信息
            parts = path.split("/")
            if len(parts) < 3:
                error_files.append(f"{path} - 路径格式错误")
                continue
                
            bioproject = parts[-4]
            sample_id = parts[-2]
            colname = f"{bioproject}_{sample_id}"

            # 读取文件，跳过前2行
            # 使用awk命令直接提取第5列，更高效
            import subprocess
            result = subprocess.run(
                f"awk 'NR>=3 {{print $1 \"\t\" $5}}' '{path}'",
                shell=True, 
                capture_output=True, 
                text=True
            )
            
            if result.returncode != 0:
                error_files.append(f"{path} - 读取失败: {result.stderr}")
                continue
                
            # 将结果转换为DataFrame
            from io import StringIO
            df_temp = pd.read_csv(StringIO(result.stdout), sep="\t", header=None, 
                                names=["bin", colname])
            
            # 设置索引
            df_temp = df_temp.set_index("bin")
            
            # 合并到总表
            if df_final is None:
                df_final = df_temp
            else:
                df_final = df_final.join(df_temp, how="outer")
                
            success_count += 1
            
        except Exception as e:
            error_files.append(f"{path} - 错误: {str(e)}")
            continue

# 处理完成后的操作
if df_final is not None:
    # 填充NaN值为0
    df_final = df_final.fillna(0)
    
    # 输出合并文件
    df_final.to_csv(output_file, sep="\t")
    
    print(f"\n合并完成!")
    print(f"成功处理: {success_count} 个文件")
    print(f"输出文件: {output_file}")
    print(f"矩阵大小: {df_final.shape[0]} 行(bins) x {df_final.shape[1]} 列(样本)")
else:
    print("没有成功处理任何文件")

# 记录错误文件
if error_files:
    print(f"处理失败: {len(error_files)} 个文件")
    with open(error_log, "w") as f:
        for error in error_files:
            f.write(error + "\n")
    print(f"错误日志: {error_log}")

print("程序执行完毕!")
