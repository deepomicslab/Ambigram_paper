## Introduction
[![DOI](https://zenodo.org/badge/428970777.svg)](https://zenodo.org/badge/latestdoi/428970777)

This is a repository for experiments on Ambigram. Here, we aim to benchmark the efficacy of Ambigram to decipher various BFB events on Illumina pair-end (PE) reads, Oxford Nanopore (ONT) long reads, Pacific Biosciences (PB) long reads, 10x Genomics linked-reads with varying tumor purity and sequencing depth. 

## Steps of in silico experiments

1. Benchmark based on simulated data: we simulated 4 sets of data to test the efficacy of Ambigram to decipher different BFB paths, including fold-back inversions, deletions and translocations. 
1. Construct a BFB path from a set of test data (sv and seg files) with Ambigram.
1. Generate a BFB fasta file (base sequence) with reference to hg38.fa. 
1. Simulate sequencing reads on the BFB fasta with different simulators for PE, PB, ONT, and 10x, respectively.
1. Align the simulated reads with the Homo sapiens (human) genome reference (e.g., hg38.fa).
1. Extract SV (structural variant) information from the BAM file derived from alignment.
1. Reconstruct the BFB path with Ambigram and new SV information. 

**Note: there is a .sh file in each directory for the specific experiment. **

## Steps of reconstructing BFB path using curated SV

1. Extract partial reads from the whole genome bam file.
2. (optional) Test the effect of tumor purity and sequencing depth.
   1. Merge the normal sample with the tumor (with SVs) sample in a ratio.
   1. Subsample the clipped bam file to generate bam files with different depths.
3. Call SVs from the clipped bam file with some tools, e.g., svaba or sniffles.
4. Convert the vcf file into a sv file.
5. Generate a lh file with the sv file and seg file.
6. Reconstruct the BFB path with Ambigram and new SV information. 

## Prerequisites

- SVAS ([https://github.com/paprikachan/SVAS](https://github.com/paprikachan/SVAS))
- localHap (to be released)
- wgsim ([https://github.com/lh3/wgsim](https://github.com/lh3/wgsim))
- bwa ([https://github.com/lh3/bwa](https://github.com/lh3/bwa))
- samtools ([https://github.com/samtools/samtools](https://github.com/samtools/samtools))
- svaba ([https://github.com/walaj/svaba](https://github.com/walaj/svaba))
- pbsim2 ([https://github.com/yukiteruono/pbsim2](https://github.com/yukiteruono/pbsim2))
- ngmlr ([https://github.com/philres/ngmlr](https://github.com/philres/ngmlr))
- sniffles ([https://github.com/fritzsedlazeck/Sniffles](https://github.com/fritzsedlazeck/Sniffles))
- LRSIM ([https://github.com/aquaskyline/LRSIM](https://github.com/aquaskyline/LRSIM))
