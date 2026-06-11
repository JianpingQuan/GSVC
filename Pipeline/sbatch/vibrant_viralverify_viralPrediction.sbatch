#!/bin/bash
#SBATCH --job-name=vb_vv
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=240G
#SBATCH --partition=amd_512q

source ~/.bashrc
conda activate vibrant
export PYTHONUNBUFFERED=1


echo "[INFO] Running viral prediction for $sample"

mkdir -p "$logdir/$bioproject/vibrant"
LOG_OUT="$logdir/$bioproject/vibrant/${sample}.gene.out"
LOG_ERR="$logdir/$bioproject/vibrant/${sample}.gene.err"

bash /public1/home/m6l600415/Public/pipeline/modules/vibrant_viralverify_viralPrediction.sh "$sample" "$bioproject" > "$LOG_OUT" 2> "$LOG_ERR"
