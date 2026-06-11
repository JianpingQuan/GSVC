#!/bin/bash

sample=$1
bioproject=$2

source /work/home/acq6mmxy68/pipeline/config/config.param

bindir="${OUTDIR}/${bioproject}/binning/${sample}"
checkm_outdir="${OUTDIR}/${bioproject}/checkm/${sample}"
checkm_db="/work/home/acq6mmxy68/miniconda3/envs/metagenome_env/checkm_data"
export CHECKM_DATA_PATH="$checkm_db"

mkdir -p "$checkm_outdir"

# 执行 CheckM
checkm lineage_wf -t 20 -x fa "$bindir" "$checkm_outdir"

