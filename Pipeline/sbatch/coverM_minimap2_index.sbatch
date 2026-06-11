#!/bin/bash
#SBATCH --job-name=mm2_index
#SBATCH --output=%j.out         # 将日志输出到以作业 ID 命名的文件，不要用 /dev/null
#SBATCH --error=%j.err          # 将错误输出到以作业 ID 命名的文件
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=36      # 匹配你申请的 36 核
#SBATCH --mem=90G               # 90G 内存足够构建 30G 序列的索引
#SBATCH --partition=xahcnormal

#------------------- 环境准备 -------------------
source ~/.bashrc
conda activate metagenome_env

# 定义路径变量（修改为你实际的绝对路径）
INDEX_DIR="/work/home/acq6mmxy68/dRep_strain_final_genomes/Minimap2_Index"
FASTA_IN="/work/home/acq6mmxy68/dRep_strain_final_genomes/processed_MAGs/all_MAGs_combined.fasta"
INDEX_OUT="${INDEX_DIR}/all_MAGs_combined.mmi"

# 确保索引目录存在
mkdir -p "$INDEX_DIR"

#------------------- 执行索引构建 -------------------
echo "开始构建索引: $(date)"

# 参数优化说明：
# -t 36: 使用全部 36 个线程加快 k-mer 计算
# -I 8G: [重要] 将索引分块大小设为 8G（默认 4G）。
#        因为你申请的内存是 90G，增大这个值可以显著减少后续比对时的磁盘换入换出，
#        从而极大提升 CoverM 计算丰度的效率。
# -d: 指定输出索引文件

minimap2 -t 36 -I 8G -d "$INDEX_OUT" "$FASTA_IN"

echo "索引构建完成: $(date)"
