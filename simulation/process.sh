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

# example (find a bfb path according to the sample data)
generate_lh test1_sv.txt test1_seg.txt test1 0
find_BFB test1.1.lh test1
