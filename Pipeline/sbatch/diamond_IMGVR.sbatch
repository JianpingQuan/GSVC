#!/bin/bash
#SBATCH --job-name=diamond
#SBATCH --output=/public1/home/m6l600415/faa_query_split/slurm_%j.out
#SBATCH --error=/public1/home/m6l600415/faa_query_split/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --partition=amd_512q

# 激活环境的稳健方式
source ~/.bashrc
conda activate base

# 取消 Python 缓冲输出
export PYTHONUNBUFFERED=1

diamond blastp \
    --query /public1/home/m6l600415/faa_query_split/query_split.${i}.faa \
    --db /public1/home/m6l600415/IMGVR_all_proteins-high_confidence.dmnd \
    --out /public1/home/m6l600415/faa_query_split/query_split.${i}_IMG_VR.tsv \
    --outfmt 6 qseqid sseqid pident length qlen evalue bitscore \
    --threads 32 \
    --max-target-seqs 1 \
    --evalue 1e-5

