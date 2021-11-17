#!/bin/bash

# extract regional data from a whole genome bam file
# para:
# 1: whole genome bam file
# 2: region file (in bed format)
# 3: prefix
function extraction()
{
    samtools view -H ${1} > ${3}.sam
    while read line;do
        samtools view ${1} $line >> ${3}.sam
    done < $2
    samtools sort -O BAM ${3}.sam -o ${3}.bam --threads 8
    samtools index ${3}.bam
    rm -f ${3}.sam
}

# example
extraction COLO829T.PE.Bam ../regions.txt PE
# Then, use other simulation scripts like PE_process.sh to do the experiment.