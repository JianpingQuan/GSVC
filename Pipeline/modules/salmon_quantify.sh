#!/bin/bash

# 参数
sample=$1
bioproject=$2

# 加载配置参数
INDEX="/work/home/acq6mmxy68/resource/salmonIndex"
Reads="/work/home/acq6mmxy68/Public/results/${bioproject}/nohost/${sample}"
THREADS=64

salmon quant -i $INDEX \
    -l A \
    -p $THREADS \
    --meta \
    -1 $Reads/${sample}_unmap_1.fq.gz \
    -2 $Reads/${sample}_unmap_2.fq.gz \
    -o /work/home/acq6mmxy68/Public/results/${bioproject}/salmon/${sample}.quant
