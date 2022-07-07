#! /usr/bin/env bash

#BSUB -J minimap2 
#BSUB -o logs/mm2.%J.out
#BSUB -e logs/mm2.%J.err
#BSUB -n 16
#BSUB -R "select[mem>100] rusage[mem=100] span[hosts=1]"
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

refn="RNA" #because we align to several references, helpful to add name in outputs
ref="~/ref/gffread_transcripts_from_v5_transcriptome_13AUG20218.fa"

# align the reads to the reference genome in a splice aware, forward strand only fashion
# this assumes you've catted all the fastqs from the run together (with any desired QC filtering) into an allreads.fastq file
# for our purposes, we've opted not to pre-filter based on read quality

minimap2 -ax splice -uf \
    -k14 ${ref} \
    allreads_combined.fastq | samtools sort -o mm2.${refn}.bam 

# generate normalization factor for genomecov to normalize to 5' end counts per million reads 
tmpscale=$(samtools view -c mm2.${refn}.bam | awk '{print 1000000/$1}')

# map 5' ends for + and - strands
bedtools genomecov -ibam mm2.${refn}.bam -5 -strand - -scale ${tmpscale} -bg > 5p_neg_${refn}.bg
bedtools genomecov -ibam mm2.${refn}.bam -5 -strand + -scale ${tmpscale} -bg > 5p_pos_${refn}.bg

# make bigwigs
bedGraphToBigWig 5p_pos_${refn}.sorted.bg ${chromsize} 5p_pos_${refn}.bw
bedGraphToBigWig 5p_neg_${refn}.sorted.bg ${chromsize} 5p_neg_${refn}.bw

