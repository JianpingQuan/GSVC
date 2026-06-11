#!/bin/bash
#SBATCH --job-name=prodigal
#SBATCH --output=/public1/home/m6l600415/resource/slurm_%j.out
#SBATCH --error=/public1/home/m6l600415/resource/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --partition=amd_512q

# 激活环境的稳健方式
ource ~/.bashrc
conda activate vibrant

# 取消 Python 缓冲输出
export PYTHONUNBUFFERED=1

prodigal -i /public1/home/m6l600415/resource/high_medium_quality_viral.Final.fasta \
         -a /public1/home/m6l600415/resource/high_medium_quality_viral.Final.proteins.faa \
         -d /public1/home/m6l600415/resource/high_medium_quality_viral.Final.genes.ffn \
         -o /public1/home/m6l600415/resource/high_medium_quality_viral.Final.prodigal_output.txt \
         -p meta

