#!/bin/bash

asm=$1
fastq1=$2
fastq2=$3
path=$4
threads=$5

if [ ! -d $path ]; then
    mkdir -p $path
fi

# map back reads to assembly to get coverage
bowtie2-build $asm $asm > $path/bowtie_build.out
bowtie2 -p $threads -x $asm -1 $fastq1 -2 $fastq2 2> $path/bowtie.out | \
samtools view -@ $threads -O BAM | \
samtools sort -@ $threads -O BAM -o $asm.bam

# calculate average coverage and generate cvg file
if [ -f $asm.bam ]; then
    coverage_eval.py $asm.bam $asm $path
    if [ $? -ne 0 ]; then
        echo " coverage eval failed with code " $?
        exit 1
    fi
else
    echo "bowtie mapping failed, $asm.bam missing"
    exit 1
fi


assembly-stats -u  $asm | cut -f11 > $path/score
# ALE score, currently has conflicts in bioconda
#ALE $asm.sorted.bam $asm $path/ALE.out > $path/ALE.log
#head -n 1 $path/ALE.out | cut -f3 -d ' ' >> $path/score






