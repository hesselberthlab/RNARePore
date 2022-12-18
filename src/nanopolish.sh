#BSUB -J nanopolish 
#BSUB -o logs/nanopolish.%J.out
#BSUB -e logs/nanopolish.%J.err
#BSUB -n 1
#BSUB -R "select[mem>100] rusage[mem=100] span[hosts=1]"
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# Sample Nanopolish script for generating intermediate files on cluster
# Here, we're resquiggling a short sequence aligned with BWA (which barfs at Us, had to convert to Ts)

# set conda environment
#conda activate nanopolish

# fix path to plugin
export HDF5_PLUGIN_PATH=/beevol/home/whitel/src/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

# build an index
nanopolish index -d fast5_pass_singles -s out/sequencing_summary.txt allreads.fastq

# combine signal and sequence data
nanopolish eventalign \
    --reads allreads.fastq \
    --bam allreadsT.bwa.sc1splint.bam \
    --genome splint.fa \
    --scale-events > eventalign_output.txt
