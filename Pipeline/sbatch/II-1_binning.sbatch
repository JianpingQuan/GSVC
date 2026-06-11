#!/bin/bash
#SBATCH --job-name=binning
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=xahcnormal
#SBATCH --export=ALL

#-------------------环境准备-------------------
source /work/home/acq6mmxy68/software/metabat/install/env.sh

# ------------------ 参数处理 ------------------
mkdir -p "$logdir/$bioproject/binning"
LOG_OUT="$logdir/$bioproject/binning/${sample}.binning.out"
LOG_ERR="$logdir/$bioproject/binning/${sample}.binning.err"

# ------------------ 执行任务 ------------------
{
  echo "[INFO] Running binning for sample: $sample"
  bash /work/home/acq6mmxy68/pipeline/modules/II-1_binning.sh "$sample" "$bioproject"
  echo "[INFO] Finished binning for $sample"
} > "$LOG_OUT" 2> "$LOG_ERR"
