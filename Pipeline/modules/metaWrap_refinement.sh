#!/bin/bash

sample=$1
bioproject=$2

#
bin_dir="/work/home/acq6mmxy68/Public/results/${bioproject}/binningWrap/${sample}/bin"
output_dir="/work/home/acq6mmxy68/Public/results/${bioproject}/binningWrap/${sample}/bin/refinement"


# 创建输出目录
mkdir -p "$output_dir"

# 执行 metaWrap 进行 binning
metawrap bin_refinement \
  -o "$output_dir" \
  -t 32 -c 50 -x 10 \
  -A "$bin_dir/concoct_bins" \
  -B "$bin_dir/maxbin2_bins" \
  -C "$bin_dir/metabat2_bins"
