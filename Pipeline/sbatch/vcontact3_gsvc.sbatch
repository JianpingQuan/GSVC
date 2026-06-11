#!/bin/bash
#SBATCH --job-name=vc3_gsvc
# 勿用 /dev/null：作业在脚本第一条之前就失败时，Slurm 仍能写 %j 日志。
# 目录必须在「提交前」已存在：login 上先执行 mkdir -p /work/home/acq6mmxy68
# （直接写到家目录下，避免 sbatch_logs 未建导致作业无法启动）
#SBATCH --output=/work/home/acq6mmxy68/vc3_slurm_%j.out
#SBATCH --error=/work/home/acq6mmxy68/vc3_slurm_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --partition=xahcnormal

# vConTACT3：GSVC strict00 PC-bucket 代表集
# 用户全流程日志（conda + vcontact3）：${OUT_ROOT}/${JOB_TAG}.${SLURM_JOB_ID}.log
# Slurm 自带日志：${WORKDIR}/vc3_slurm_${SLURM_JOB_ID}.out / .err

#------------------- 用户参数（按需修改）-------------------
WORKDIR="/work/home/acq6mmxy68"
OUT_ROOT="${WORKDIR}"
JOB_TAG="vcontact3_run2"

PROTEINS="${WORKDIR}/GSVC_strict00_PCbucket_top1.faa"
G2G="${WORKDIR}/GSVC_strict00_PCbucket_top1_gene2genome.parquet"
NUCL_LEN="${WORKDIR}/GSVC_strict00_PCbucket_top1_genome_lengths.parquet"
DB_PATH="${WORKDIR}/resource/vcontact3"
DB_VERSION="230"
OUTPUT_DIR="${WORKDIR}/${JOB_TAG}"

THREADS="${SLURM_CPUS_PER_TASK:-32}"

#------------------- 尽早落盘：任何后续失败都能看见 -------------------
# 注意：必须在 conda / set -e 之前执行，否则 conda 失败时永远不会执行到原来的 { } > log
mkdir -p "${OUT_ROOT}"
LOG_ALL="${OUT_ROOT}/${JOB_TAG}.${SLURM_JOB_ID:-nojid}.log"
exec >> "${LOG_ALL}" 2>&1
echo "========== $(date -Is) job start =========="
echo "[INFO] SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "[INFO] host=$(hostname) pwd=$(pwd)"
echo "[INFO] bash=$(command -v bash) BASH_VERSION=${BASH_VERSION:-}"

# 勿用 set -u：conda 的 deactivate.d/*.sh 会引用未定义的 CONDA_BACKUP_*，与 nounset 冲突导致激活失败
set -eo pipefail
if [ "${DEBUG_VC3_SBATCH:-0}" = 1 ]; then set -x; fi

#------------------- Conda（与登录节点路径一致时可写死一行 source）-------------------
_conda_sh=""
for _p in \
  "${HOME}/miniconda3/etc/profile.d/conda.sh" \
  "${WORKDIR}/miniconda3/etc/profile.d/conda.sh" \
  "${HOME}/anaconda3/etc/profile.d/conda.sh" \
  "${HOME}/.conda/etc/profile.d/conda.sh"; do
  if [ -f "${_p}" ]; then
    _conda_sh="${_p}"
    break
  fi
done
if [ -z "${_conda_sh}" ]; then
  echo "[ERROR] 未找到 conda.sh。请在登录节点执行: ls ~/miniconda3/etc/profile.d/conda.sh"
  echo "[ERROR] 若 conda 在其它路径，把本脚本里 for 循环改成你的 conda.sh 绝对路径。"
  exit 1
fi
# conda activate 会跑 deactivate.d；在 nounset 下易触发未定义变量错误，故全程不用 set -u
set +u
# shellcheck source=/dev/null
source "${_conda_sh}"
conda activate vcontact3

command -v vcontact3
vcontact3 version

#------------------- 无图形界面（Slurm 计算节点无 X11）-------------------
# 导出 newick 等时 ETE3/Qt 默认连 xcb 会 Abort；强制 Qt 离屏 + Matplotlib 非交互
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"
export MPLBACKEND="${MPLBACKEND:-Agg}"
unset DISPLAY

#------------------- 执行任务 ------------------
echo "[INFO] cpus: ${THREADS}"
echo "[INFO] proteins: ${PROTEINS}"
echo "[INFO] gene2genome: ${G2G}"
echo "[INFO] len-nucleotide: ${NUCL_LEN}"
echo "[INFO] db-path: ${DB_PATH} db-version: ${DB_VERSION}"
echo "[INFO] output: ${OUTPUT_DIR}"

vcontact3 run \
  --proteins "${PROTEINS}" \
  --gene2genome "${G2G}" \
  --len-nucleotide "${NUCL_LEN}" \
  --db-domain prokaryotes \
  --db-path "${DB_PATH}" \
  --db-version "${DB_VERSION}" \
  --exports cytoscape graphml profiles newick \
  --target-rank order family \
  --target-members 3 \
  --threads "${THREADS}" \
  --reduce-memory \
  --no-progress \
  --output "${OUTPUT_DIR}" \
  --force-overwrite
echo "========== $(date -Is) job end OK =========="
