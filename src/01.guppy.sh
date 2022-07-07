#! /usr/bin/env bash

#BSUB -J guppy 
#BSUB -o logs/guppy.%J.out
#BSUB -e logs/guppy.%J.err
#BSUB -n 24
#BSUB -q rna
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# rebasecall RNA in high accuracy mode 
guppy_basecaller -i fast5_all -s out -c rna_r9.4.1_70bps_hac.cfg \
    --fast5_out --num_callers 6
