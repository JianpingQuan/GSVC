#!/bin/bash
#SBATCH --job-name=drep
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --partition=xahcnormal

source ~/.bashrc
conda activate metagenome_env

set -e

echo "[INFO] Job started at $(date) on $(hostname)"

# ========================
# 路径配置（你可以根据需要修改）
# ========================
DREP_OUT="/work/home/acq6mmxy68/Public/results/"$bioproject"/drep/drep_out"
GENOME_DIR="/work/home/acq6mmxy68/Public/results/"$bioproject"/drep/drep_input"
THREADS=24

# ========================
# 创建输出目录与日志目录
# ========================
mkdir -p "$DREP_OUT"

# ========================
# 运行 dRep
# ========================
echo "[INFO] Running dRep dereplicate..."

dRep dereplicate "$DREP_OUT" \
  -g ${GENOME_DIR}/*.fa \
  -p "$THREADS" \
  -comp 95 -con 5 \
  -sa 0.99 -nc 0.3

echo "[✅] dRep dereplicate completed at $(date)"
