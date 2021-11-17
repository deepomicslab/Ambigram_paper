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

# simulate 10x linked-reads sequencing
# para:
# 1: coverage (30)
# 2: read_len (150)
# 3: insertion_size (350)
# 4: error_rate 
# 5: prefix
# 6: input fasta
# 7: molecule_len (50)
function tenx_simulation()
{
    fasta_length=$(cat $6 | grep -v '^>' | wc -lc | awk '{print $2 - $1}')
    let "10X_reads = $fasta_length*$1/300"
    10X_reads=$(echo "$10X_reads 1000000" | awk '{printf "%.6f", $1/$2}')
    let "10X_length = $fasta_length*$1/2"
    partition_lower=$(echo "$10X_length $7" | awk '{printf "%.6f", $1/$2/1000}')
    let "partition_lower = $(printf %.0f $partition_lower)+1"
    partition_lower=$(echo "$partition_lower 1000" | awk '{printf "%.3f", $1/$2}')

    # note: for tumor (with SVs) ref: -r  for normal ref: -g
    ~/LRSIM/simulateLinkedReads.pl -r $6 -p $5 -d 1 -1 $fasta_length  -i $3 -x $10X_reads -f $7 -t $partition_lower -m 1 -o -4 1 -7 1 -e $4
}

# alignment for 10x linked-reads 
# param:
# 1: sample name (prefix)
# 2: fastq directory (default: ./)
# 3: reference directory
function tenx_alignment()
{
    longranger align \
        --id=$1 \
        --reference=$3 \
        --fastqs=$2 \
        --sample=$1
}

# SV calling (generate .vcf file)
# param:
# 1: reference fasta
# 2: bam file
# 3: prefix 
function call_sv()
{
    svaba run -p 4 -t $2 -a $3 -G $1
}

# generate sv file based on svaba .vcf
# 1: vcf file
# 2: output prefix
function get_sv()
{
    python ~/SVAS/scripts/parse_svaba.py call --sv_fn=$1
    python ~/localhap/localhaptgs/script/main.py vcf2sv -i ${1}.txt -o ${2}_sv.txt
}

# example
call_sv hg19.fa 10X.bam 10X
get_sv 10X.unfiltered.sv.vcf 10X
generate_lh 10X_sv.txt ../COLO829_1_seg.txt sample1_10X 1
find_BFB sample1_10X.1.lh sample1