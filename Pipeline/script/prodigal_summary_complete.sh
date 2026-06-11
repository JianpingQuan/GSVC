#!/bin/bash

out="prodigal_summary_complete.tsv"
echo -e "Sample\tNum_seqs\tSum_len\tMin_len\tAvg_len\tMax_len" > $out

for f in /work/home/acq6mmxy68/pipeline/logs/*/prodigal/*.gene.out; do
    sample=$(basename "$f" .gene.out)   # 提取样本名
    stats=$(grep -B1 'completed' "$f" | awk 'NR==1 {OFS="\t"; print $4, $5, $6, $7, $8}'|sed 's/,//g')
    echo -e "${sample}\t${stats}" >> $out
done
