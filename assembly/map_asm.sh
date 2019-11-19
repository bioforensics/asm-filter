#!/bin/bash
#
# script takes assembly and maps reads to assembly
#

#output dir
dir=$1
map_ctg=$2
read1=$3
read2=$4

prefix=$(basename $map_ctg .fasta)

if [ ! -d $dir ]; then
    mkdir -p $dir
fi


if [ -f $map_ctg ]; then

    bowtie2-build $map_ctg $map_ctg > $dir/bowtie_build.out
    
    bowtie2 -p 32 -x $map_ctg -1 $read1 -2 $read2 2> $dir/bowtie2.log | \
    samtools view -@ 32 -O BAM | samtools sort -@ 32 -O BAM -o $dir/${prefix}.bam
    
    # Extract aligned reads
    samtools view -@ 32 -F4 -o $dir/${prefix}.aligned_reads.bam -O BAM $dir/${prefix}.bam
    #$samtools view -f4 -o $dir.${prefix}.bowtie2.sorted.unaligned.bam -O BAM $dir.${prefix}.bowtie2.sorted.bam

    # Extract both pairs unaligned reads
    samtools view -@ 32 -f12 -F 256 -O BAM -o $dir/${prefix}.both_unaligned.bam $dir/${prefix}.bam

    
    # mapped and unmapped read extraction
    bedtools bamtofastq -i $dir/${prefix}.aligned_reads.bam -fq $dir/${prefix}.aligned_reads.fq
    
    bedtools bamtofastq -i $dir/${prefix}.both_unaligned.bam -fq $dir/${prefix}.both_unaligned.R1.fq -fq2 $dir/${prefix}.both_unaligned.R2.fq
    
    # convert aligned reads to fasta
    seqtk seq -A $dir/${prefix}.aligned_reads.fq > $dir/${prefix}.aligned_reads.fa
    
    fastaLengths.pl $dir/${prefix}.aligned_reads.fa | cut -f 1 | stats.pl - > $dir/${prefix}.aligned_reads.stats
    fastaLengths.pl $map_ctg > $dir/mapping.contig.fasta.len
    
else
    echo "$map_ctg doesn't exist!"
    exit
fi

