#!/bin/bash
#SBATCH --job-name=iphop_part
#SBATCH --output=/public1/home/m6l600415/resource/slurm_%j.out
#SBATCH --error=/public1/home/m6l600415/resource/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --partition=amd_512q

# 激活环境的稳健方式
source /public1/home/m6l600415/miniforge3/bin/activate iphop38
source /public1/soft/modules/module.sh
module load gcc/7.3.0-wzm
ulimit -n 655350

# 取消 Python 缓冲输出
export PYTHONUNBUFFERED=1

# 注意：去掉末尾的 &，让 SLURM 正确监控进程
# 变量 i 已经在 sbatch --export 中传递进来
iphop predict --fa_file /public1/home/m6l600415/resource/part_${i}.fasta \
              --out_dir /public1/home/m6l600415/resource/out_part_${i} \
              --db_dir /public1/home/m6l600415/software/iphop/Jun_2025_pub_rw/ \
              --num_threads 32 > /public1/home/m6l600415/resource/log_part_${i}.txt 2>&1

# 检查退出状态
if [ $? -eq 0 ]; then
    echo "Part ${i} finished successfully"
else
    echo "Part ${i} failed" >&2
    exit 1
fi
