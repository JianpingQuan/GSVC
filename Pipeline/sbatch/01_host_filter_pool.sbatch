#!/bin/bash
#SBATCH --job-name=host_fliter
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=28G
#SBATCH --partition=xahcnormal


#-------------------环境准备-------------------
source ~/.bashrc
conda activate metagenome_env


# ------------------ 参数处理 ------------------
echo "[INFO] Loading config from: $contig"
source "$contig" || { echo "Failed to source $contig"; exit 1; }

mkdir -p "$logdir/$bioproject/hostFilter"
LOG_OUT="$logdir/$bioproject/hostFilter/${sample}.hostFilter.out"
LOG_ERR="$logdir/$bioproject/hostFilter/${sample}.hostFilter.err"


# ------------------ 执行任务 ------------------
{
  echo "[INFO] Running hostFilter for sample: $sample"
  echo "[INFO] fq1: $fq1"
  echo "[INFO] fq2: $fq2"

  bash /work/home/acq6mmxy68/pipeline/modules/01_host_filter_pool.sh $sample $bioproject $fq1 $fq2
  
  echo "[INFO] Finished host filtering for $sample"
} > "$LOG_OUT" 2> "$LOG_ERR"
