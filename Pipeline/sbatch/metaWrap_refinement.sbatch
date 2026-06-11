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
source ~/.bashrc
conda activate metawrap

# ------------------ 参数处理 ------------------
mkdir -p "$logdir/$bioproject/binningWrap"
LOG_OUT="$logdir/$bioproject/binningWrap/${sample}.refinement.out"
LOG_ERR="$logdir/$bioproject/binningWrap/${sample}.refinement.err"

# ------------------ 执行任务 ------------------
{
  echo "[INFO] Running binning for sample: $sample"

  bash /work/home/acq6mmxy68/pipeline/modules/metaWrap_refinement.sh "$sample" "$bioproject"
 
 echo "[INFO] Finished binning for $sample"
} > "$LOG_OUT" 2> "$LOG_ERR"
