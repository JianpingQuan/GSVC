#!/bin/bash

#BATCH --job-name=coverM
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --partition=xahcnormal
#SBATCH --exclusive

#-------------------环境准备-------------------
source ~/.bashrc
conda activate metagenome_env 

# ------------------ 参数处理 ------------------
mkdir -p "$logdir/$bioproject/viral_coverm"
LOG_OUT="$logdir/$bioproject/viral_coverm/${sample}.coverm.database1.out"
LOG_ERR="$logdir/$bioproject/viral_coverm/${sample}.coverm.database1.err"

# ------------------ 执行任务 ------------------
{
  DATE=$(date "+%Y-%m-%d %H:%M:%S")

  echo "[INFO] Running CoverM for sample: $sample - $DATE"

  bash /work/home/acq6mmxy68/pipeline/script/viral_coverM.sh "$sample" "$bioproject"

  DATE=$(date "+%Y-%m-%d %H:%M:%S")

  echo "[INFO] Finished CoverM for $sample - $DATE"
} > "$LOG_OUT" 2> "$LOG_ERR"
