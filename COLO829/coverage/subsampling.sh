#!/bin/bash

# This script shows how to subsample a bam file to simulate sequencing reads with different depths (coverage).
# For the specific experiments, please refer to the scripts of different simulators.

# subsampling for bam file
# 1: subsampling ratio
# 2: input bam file
# 3: prefix
function subsampling()
{
    samtools view -s $1 -b $2 > ${3}.sam
    samtools sort -O BAM ${3}.sam -o ${3}.bam --threads 8
    samtools index ${3}.bam
}

# example (subsample 20x depth from a bam file with 30x depth)
subsampling 0.67 ./PE.bam PE_20x