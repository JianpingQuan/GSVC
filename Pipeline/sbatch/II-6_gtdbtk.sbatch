#!/bin/bash
#SBATCH --job-name=gtdbtk_classify
#SBATCH --error=/work/home/acq6mmxy68/Public/results/00_CombinedAnalysis/gtdbtk_out/gtdbtk_classify.%j.err
#SBATCH --output=/work/home/acq6mmxy68/Public/results/00_CombinedAnalysis/gtdbtk_out/gtdbtk_classify.%j.out
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH -N 1
#SBATCH --partition=xahcnormal
#SBATCH --exclusive
#SBATCH --constraint="64core"

source ~/.bashrc
conda activate gtdbtk_py310 

set -e

echo "[INFO] Job started at $(date) on $(hostname)"

# ========================
# 路径配置（你可以根据需要修改）
# ========================
Input_Dir="/work/home/acq6mmxy68/Public/results/00_CombinedAnalysis/drep_out2/dereplicated_genomes/${group}"
Output_Dir="/work/home/acq6mmxy68/Public/results/00_CombinedAnalysis/gtdbtk_out/${group}"
THREADS=64

# ========================
# 创建输出目录与日志目录
# ========================
mkdir -p "$Output_Dir"

# ========================
# 运行 gtdbtk
# ========================
echo "[INFO] Running gtdbtk classify_wf..."

gtdbtk classify_wf --genome_dir "${Input_Dir}" --out_dir "${Output_Dir}" --extension fa --cpus $THREADS --skip_ani_screen  --pplacer_cpus 1 --scratch_dir /work/home/acq6mmxy68/Public/results/00_CombinedAnalysis/gtdbtk_out/

echo "[INFO] gtdbtk classify_wf completed at $(date)"
