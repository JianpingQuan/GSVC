#!/bin/bash
#SBATCH --job-name=virsorter
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=xahcnormal
#SBATCH --export=ALL

# ------------------ Conda环境激活 ------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate vs2

# ------------------ 参数处理 ------------------
mkdir -p "$logdir/$bioproject/virsorter"
LOG_OUT="$logdir/$bioproject/virsorter/${sample}.virsorter.out"
LOG_ERR="$logdir/$bioproject/virsorter/${sample}.virsorter.err"

# ------------------ 执行任务 ------------------
{
  echo "[INFO] Running virsorters for sample: $sample at: $(date)"
  bash /work/home/acq6mmxy68/pipeline/modules/virsorter.sh "$sample" "$bioproject"
  echo "[INFO] Finished virsorter for sample: $sample at: $(date)"
} > "$LOG_OUT" 2> "$LOG_ERR"
