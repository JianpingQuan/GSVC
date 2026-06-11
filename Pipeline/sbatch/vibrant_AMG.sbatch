#!/bin/bash
#SBATCH --job-name=vb_vv
#SBATCH --output=/public1/home/m6l600415/Public/amg/vibrant_%j.out
#SBATCH --error=/public1/home/m6l600415/Public/amg/vibrant_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --partition=amd_512q

source ~/.bashrc
conda activate vibrant

ulimit -n 655350
export PYTHONUNBUFFERED=1

#INPUT_FASTA="/public1/home/m6l600415/resource/high_quality_viral_rep_seq.fasta"
INPUT_FASTA="/public1/home/m6l600415/resource/medium_quality_vOTU.fasta"
OUT_FOLDER="/public1/home/m6l600415/Public/amg"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input file $INPUT_FASTA not found!"
    exit 1
fi

VIBRANT_run.py \
    -i "$INPUT_FASTA" \
    -folder "$OUT_FOLDER" \
    -t 32 \
    -no_plot

# 检查退出状态
if [ $? -eq 0 ]; then
    echo "VIBRANT analysis finished successfully at $(date)"
else
    echo "VIBRANT analysis failed at $(date)" >&2
    exit 1
fi
