#!/bin/bash

# This script shows how to merge the normal sample with the tumor (SVs) sample to simulate samples with different tumor purity.
# For the specific experiments, please refer to the scripts of different simulators.

# generate a sample with a specified tumor purity
# 1: normal BAM
# 2: normal purity
# 3: tumor BAM
# 4: tumor purity
# 5: prefix
function merge()
{
    samtools view -b -s $2 $1 > temp1.bam
    samtools view -b -s $4 $3 > temp2.bam
    samtools merge ${5}.bam temp1.bam temp2.bam -f
    samtools index ${5}.bam
    rm -f temp1.bam
    rm -f temp2.bam
}

# example (generate a bam file with 75% tumor purity)
merge normal.bam 0.25 PE.bam 0.75 PE_75x
