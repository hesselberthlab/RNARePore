#! /usr/bin/env bash

#BSUB -J ont-preprocess
#BSUB -o logs/ont-preprocess.%J.out
#BSUB -e logs/ont-preprocess.%J.err
#BSUB -n 16
#BSUB -R "span[hosts=1]"

source activate tombo

set -o nounset -o pipefail -o errexit -x

# locations after rebasecalling
fast5_data="/out/workspace"
fastq_data="/out/pass"
singles="fast5_pass_singles"
path="$HOME/path/to/run/data" # edit this to point to each run

export H5PY_DEFAULT_READONLY=1

multi_to_single_fast5 --input_path ${path}/$fast5_data \
  --save_path $singles \
  --threads 16

tombo preprocess annotate_raw_with_fastqs \
  --overwrite --fast5-basedir $singles \
  --processes 16 \
  --fastq-filenames ${path}/${fastq_data}/*.fastq

tombo resquiggle $singles $ref \
    --processes 16 \
    --num-most-common-errors 5 \
    --ignore-read-locks
