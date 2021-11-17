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

# simulate Paired-End sequencing
# para:
# 1: coverage (30) 
# 2: read_len (150)
# 3: insertion_size (350)
# 4: error_rate 
# 5: prefix
# 6: input fasta
function PE_simulation()
{
    # coverage = (read count * read length ) / total genome size
    # For PE sequencing, read count doubles
    fasta_length=$(cat $6 | grep -v '^>' | wc -lc | awk '{print $2 - $1}')
    let "pe_reads = $fasta_length*$1/$2"
   
    wgsim -1 $2 -2 $2 -d $3 -N $pe_reads -r 0 -R 0 -X 0 -e $4 $6 ${5}_r1.fq ${5}_r2.fq
    cat ${5}_r1.fq | gzip -c > ${5}_r1.fq.gz
    cat ${5}_r2.fq | gzip -c > ${5}_r2.fq.gz
    rm ${5}_r1.fq ${5}_r2.fq
}

# alignment for paired-end reads 
# param:
# 1: forward reads
# 2: reverse reads
# 3: prefix 
# 4: reference
function PE_alignment()
{
    bwa mem -t 8 $4 $1 $2 > ${3}.sam
    samtools sort -O BAM ${3}.sam -o ${3}.bam --threads 8
    samtools index ${3}.bam
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
generate_lh ../test1_sv.txt ../test1_seg.txt test1 0
find_BFB test1.1.lh test1
get_fasta hg38.fa test1.bed test1.fa
PE_simulation 30 150 350 0 bfb1 test1.fa
PE_alignment ./bfb1_r1.fq.gz ./bfb1_r2.fq.gz bfb1 hg38.fa
call_sv hg38.fa bfb1.bam bfb1
get_sv bfb1.unfiltered.sv.vcf bfb1
generate_lh bfb1_sv.txt ../test1_seg.txt bfb1 1
find_BFB bfb1.1.lh bfb1