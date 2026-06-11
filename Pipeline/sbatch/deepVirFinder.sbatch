#!/bin/bash
#SBATCH --job-name=deepVirFinder
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=xahcnormal
#SBATCH --export=ALL

# ------------------ Conda环境激活 ------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dvf

# ------------------ 参数处理 ------------------
mkdir -p "$logdir/$bioproject/dvf"
LOG_OUT="$logdir/$bioproject/dvf/${sample}.dvf.out"
LOG_ERR="$logdir/$bioproject/dvf/${sample}.dvf.err"

# ------------------ 执行任务 ------------------
{
  echo "[INFO] Running deepVirFinder for sample: $sample at: $(date)"
  bash /work/home/acq6mmxy68/pipeline/modules/deepVirFinder.sh "$sample" "$bioproject"
  echo "[INFO] Finished deepVirFinder for sample: $sample at: $(date)"
} > "$LOG_OUT" 2> "$LOG_ERR"
