#! /usr/bin/env bash

#BSUB -J bwa
#BSUB -o logs/bwa.%J.out
#BSUB -e logs/bwa.%J.err
#BSUB -R "select[mem>100] rusage[mem=100] span[hosts=1]"
#BSUB -n 1

bwaidx="$HOME/ref/noints-consensus-ONTadapt-sacCer3-tRNAs.fa"

set -o nounset -o pipefail -o errexit -x

# upstream of this I've done a little manual preprocessing for BWA as follows:
# if Guppy performed filtering based on quality score, combine the pass and fail fastqs
# then convert all Us to Ts for BWA alignment
# can do this with the one-liner:
# awk 'NR %4 == 2{gsub("U","T")}1' allreads.fastq > allreadsT.fastq

# some LSF job server specific language here for just 1 job
samples=(
allreadsT
)

u=${samples[$(( $LSB_JOBINDEX -1 ))]}

# align to the intronless, cytoplasmic tRNome with RNA adapters added
# bwa mem parameters chosen based on PMID: 34618430
bwa mem -W 13 -k 6 -x ont2d $bwaidx \
/beevol/home/whitel/data/2021_Nanopore/tRNAseq/tpt1d_tRNAseq/20211224_1225_MN31004_FAQ63690_1e73e2c4/out/${u}.fastq
| samtools view -F4 -hu - \
| samtools sort -o ${u}.bwa.sc3.tRNome.bam

samtools index ${u}.bwa.sc3.tRNome.bam
