#!/bin/bash

sample=$1
bioproject=$2

source /work/home/acq6mmxy68/software/metabat/install/env.sh

# 路径设置
OUTDIR="/work/home/acq6mmxy68/Public/results"
megahit_outdir="${OUTDIR}/${bioproject}/megahit/${sample}"
contig="${megahit_outdir}/${sample}.contigs.fa"
depth="${OUTDIR}/${bioproject}/depth/${sample}.depth.txt"
binning_outdir="${OUTDIR}/${bioproject}/binning/${sample}"

# 创建输出目录
mkdir -p "$binning_outdir"
mkdir -p "${OUTDIR}/${bioproject}/depth"

# 计算depth文件（bowtie2 + jgi_summarize_bam_contig_depths）

/work/home/acq6mmxy68/software/bowtie2/bowtie2-build "$contig" "${megahit_outdir}/index"
/work/home/acq6mmxy68/software/bowtie2/bowtie2 -x "${megahit_outdir}/index" -1 "${OUTDIR}/${bioproject}/nohost/${sample}/${sample}_unmap_1.fq.gz" \
          -2 "${OUTDIR}/${bioproject}/nohost/${sample}/${sample}_unmap_2.fq.gz" \
          -S "${megahit_outdir}/${sample}.sam" -p 16

/work/home/acq6mmxy68/software/samtools-1.22.1/bin/samtools view -bS "${megahit_outdir}/${sample}.sam" > "${megahit_outdir}/${sample}.bam"
/work/home/acq6mmxy68/software/samtools-1.22.1/bin/samtools sort "${megahit_outdir}/${sample}.bam" -o "${megahit_outdir}/${sample}.sorted.bam"
/work/home/acq6mmxy68/software/samtools-1.22.1/bin/samtools index "${megahit_outdir}/${sample}.sorted.bam"

# 生成原始 depth 文件
depth="${OUTDIR}/${bioproject}/depth/${sample}.depth.txt"
jgi_summarize_bam_contig_depths --outputDepth "$depth" "${megahit_outdir}/${sample}.sorted.bam"


# 执行 MetaBAT2 进行 binning
metabat2 -i "$contig" -a "$depth" -o "${binning_outdir}/bin" -t 32 --minContig 1500

#若MetaBAT2成功运行，清理中间文件
if [ $? -eq 0 ]; then
  echo "MetaBAT2 succeeded. Cleaning up intermediate files..."
  rm -f "${megahit_outdir}/${sample}.sam"
  rm -f "${megahit_outdir}/${sample}.bam"
  rm -f "${megahit_outdir}/${sample}.sorted.bam"
  rm -f "${megahit_outdir}/${sample}.sorted.bam.bai"
  rm -f "${megahit_outdir}/index".*.bt2
else
  echo "MetaBAT2 failed. Intermediate files preserved for debugging." >&2
fi
