#!/bin/bash

while read line
do
    bash /work/home/acq6mmxy68/pipeline/script/02_generate_megahit_report.sh /work/home/acq6mmxy68/Public/results/"$line"/megahit
done </work/home/acq6mmxy68/pipeline/config/sampleList/bioproject.list.txt

