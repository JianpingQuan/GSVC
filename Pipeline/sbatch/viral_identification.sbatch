#!/bin/bash
#SBATCH --job-name=viralIdent
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --partition=amd_512q

set -euo pipefail
export PYTHONUNBUFFERED=1

if [[ -z "${sample:-}" || -z "${bioproject:-}" || -z "${logdir:-}" ]]; then
  echo "[ERROR] sample / bioproject / logdir not set"
  exit 1
fi

echo "[INFO] Running viral prediction"
echo "[INFO] sample=${sample}"
echo "[INFO] bioproject=${bioproject}"


mkdir -p "$logdir/$bioproject/viralIdentification"
LOG_OUT="$logdir/$bioproject/viralIdentification/${sample}.out"
LOG_ERR="$logdir/$bioproject/viralIdentification/${sample}.err"

bash /public1/home/m6l600415/Public/pipeline/modules/viralIdentification.sh "$sample" "$bioproject" > "$LOG_OUT" 2> "$LOG_ERR"
