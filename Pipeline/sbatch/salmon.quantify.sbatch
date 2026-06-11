#!/bin/bash
#SBATCH --job-name=salmon_mapping
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=160G
#SBATCH --partition=xahcnormal
##SBATCH --exclusive

source ~/.bashrc
conda activate metagenome_env

echo "==== salmon job started at $(date) ===="

mkdir -p "$logdir/$bioproject/salmon/"
LOG_OUT="$logdir/$bioproject/salmon/${sample}.salmon.out"
LOG_ERR="$logdir/$bioproject/salmon/${sample}.salmon.err"

bash /work/home/acq6mmxy68/pipeline/modules/salmon_quantify.sh "$sample" "$bioproject" > "$LOG_OUT" 2> "$LOG_ERR"

echo "==== salmon job finished at $(date) ===="
