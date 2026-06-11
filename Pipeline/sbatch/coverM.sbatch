#!/bin/bash
#SBATCH --job-name=coverM
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --partition=amd_512q

#-------------------环境准备-------------------
source ~/.bashrc
conda activate coverm_env
export PYTHONUNBUFFERED=1

# ------------------ 参数处理 ------------------
mkdir -p "$logdir/$bioproject/coverm"
LOG_OUT="$logdir/$bioproject/coverm/${sample}.coverm.out"
LOG_ERR="$logdir/$bioproject/coverm/${sample}.coverm.err"

# ------------------ 执行任务 ------------------
{
  DATE=$(date "+%Y-%m-%d %H:%M:%S")

  echo "[INFO] Running CoverM for sample: $sample - $DATE"

  bash /public1/home/m6l600415/Public/pipeline/modules/coverM.sh "$sample" "$bioproject"

  DATE=$(date "+%Y-%m-%d %H:%M:%S")

  echo "[INFO] Finished CoverM for $sample - $DATE"
} > "$LOG_OUT" 2> "$LOG_ERR"
