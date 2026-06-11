#!/bin/bash
#SBATCH --job-name=virSorter2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --partition=amd_512q

source ~/.bashrc
conda activate vs2
export PYTHONUNBUFFERED=1


echo "[INFO] Running viral prediction for $sample"

mkdir -p "$logdir/$bioproject/virsorter2"
LOG_OUT="$logdir/$bioproject/virsorter2/${sample}.gene.out"
LOG_ERR="$logdir/$bioproject/virsorter2/${sample}.gene.err"

bash /public1/home/m6l600415/Public/pipeline/modules/virsorter2_viralPrediction.sh "$sample" "$bioproject" > "$LOG_OUT" 2> "$LOG_ERR"
