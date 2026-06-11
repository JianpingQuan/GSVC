#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=xahcnormal


#-------------------环境准备-------------------
source ~/.bashrc
conda activate metagenome_env


# ------------------ 参数处理 ------------------
echo "[INFO] Loading config from: $contig"
source "$contig" || { echo "Failed to source $contig"; exit 1; }

mkdir -p "$logdir/$bioproject/qc"
LOG_OUT="$logdir/$bioproject/qc/${sample}.qc.out"
LOG_ERR="$logdir/$bioproject/qc/${sample}.qc.err"


# ------------------ 执行任务 ------------------
{
  echo "[INFO] Running QC for sample: $sample"
  echo "[INFO] fq1: $fq1"
  echo "[INFO] fq2: $fq2"

  bash /work/home/acq6mmxy68/pipeline/modules/00_raw_qc.sh "$sample" "$bioproject" "$fq1" "$fq2"

  echo "[INFO] Finished QC for $sample"
} > "$LOG_OUT" 2> "$LOG_ERR"

