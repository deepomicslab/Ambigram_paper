#!/bin/bash

# generate a .lh file (input of Onion-BFB)
# para: 
# 1: sv file
# 2: seg file
# 3: output prefix
# 4: boolean: is the sv file from simulation [0 or 1]
function generate_lh()
{
    python ~/SVAS/scripts/csv_sv.py call --sv_fn=$1 --seg_fn=$2 --sample=$3 --seg_only=1 --bfb_sv=$4
}

# run Onion-BFB
# para: 
# 1: lh file
# 2: output prefix
function find_BFB()
{
    ~/localhap/localHap --op bfb --in_lh $1 --lp_prefix $2 --verbose
}

# convert a bed file into fasta
# para: 
# 1: reference fast
# 2: bed file
# 3: output fasta
function get_fasta()
{
    python ~/localHap/localhaptgs/script/main.py bfb2fasta -i=$1 -b=$2 -t=bpseq.txt -o=${3}
}

# simulate PacBio sequencing
# para:
# 1: input fasta
# 2: prefix
# 3: depth/coverage
function PB_simulation()
{
    pbsim --prefix $2 --depth $3 --hmm_model ~/pbsim2/data/P6C4.model $1
}

# alignment for PacBio reads 
# para:
# 1: fastq file
# 2: prefix
# 3: reference fasta
function PB_alignment()
{
    ngmlr -t 16 -r $3 -q $1 -o ${2}.sam
    samtools sort -O BAM ${2}.sam -o ${2}.bam --threads 8
    samtools index ${2}.bam
}

# SV calling (generate .vcf file)
# para:
# 1: bam file
# 2: prefix 
# 3: Minimum length of SV to be reported. Default: 30bp
# 4: Minimum number of reads that support a SV to be reported. Default: 10
# 5: Maximum distance to group SV together. Sniffles estimates this parameter during runtime to group together SVs reported by different reads. Default: 1kb
function call_sv()
{
    sniffles -m $1 -v ${2}.vcf -l $3 -s $4 -d $5
}

# generate sv file based on svaba .vcf
# 1: vcf file
# 2: output prefix
function get_sv()
{
    python ~/localhap/localhaptgs/script/main.py parse_snif_vcf -i $1 -o ${2}_sv.txt
}

# example
generate_lh ../test1_sv.txt ../test1_seg.txt test1 0
find_BFB test1.1.lh test1
get_fasta hg38.fa test1.bed test1.fa
PB_simulation test1.fa bfb1 30
PB_alignment bfb1.fastq bfb1 hg38.fa
call_sv bfb1.bam bfb1 0 2 4
get_sv bfb1.vcf bfb1
generate_lh bfb1_sv.txt ../test1_seg.txt bfb1 1
find_BFB bfb1.1.lh bfb1